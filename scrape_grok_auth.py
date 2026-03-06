#!/usr/bin/env python3
"""
Authenticated Grok share URL scraper using Edge browser cookies.
Reads Edge cookies via DPAPI decryption (current user only, no external deps).
Usage: python scrape_grok_auth.py <URL>
"""

import os
import sys
import re
import json
import shutil
import sqlite3
import tempfile
import base64
import struct
import requests

import win32crypt
from Cryptodome.Cipher import AES


def get_edge_master_key():
    """Read and decrypt the Edge AES master key from Local State using DPAPI."""
    local_state_path = os.path.join(
        os.environ['LOCALAPPDATA'],
        'Microsoft', 'Edge', 'User Data', 'Local State'
    )
    with open(local_state_path, 'r', encoding='utf-8') as f:
        local_state = json.load(f)

    encrypted_key_b64 = local_state['os_crypt']['encrypted_key']
    encrypted_key = base64.b64decode(encrypted_key_b64)
    # Remove "DPAPI" prefix (first 5 bytes)
    encrypted_key = encrypted_key[5:]
    # Decrypt with DPAPI (current user)
    master_key = win32crypt.CryptUnprotectData(encrypted_key, None, None, None, 0)[1]
    return master_key


def decrypt_cookie_value(encrypted_value, master_key):
    """Decrypt a Chrome/Edge AES-GCM encrypted cookie value."""
    if encrypted_value[:3] == b'v10' or encrypted_value[:3] == b'v11':
        # AES-256-GCM: v10/v11 + 12-byte nonce + ciphertext + 16-byte tag
        nonce = encrypted_value[3:15]
        ciphertext = encrypted_value[15:]
        cipher = AES.new(master_key, AES.MODE_GCM, nonce=nonce)
        plaintext = cipher.decrypt(ciphertext[:-16])  # strip GCM tag
        return plaintext.decode('utf-8', errors='replace')
    else:
        # Legacy DPAPI-encrypted value
        try:
            return win32crypt.CryptUnprotectData(encrypted_value, None, None, None, 0)[1].decode('utf-8', errors='replace')
        except Exception:
            return ''


def get_xcom_cookies():
    """Extract auth_token and ct0 cookies for x.com from Edge's cookie database."""
    cookies_db = os.path.join(
        os.environ['LOCALAPPDATA'],
        'Microsoft', 'Edge', 'User Data', 'Default', 'Network', 'Cookies'
    )
    if not os.path.exists(cookies_db):
        raise FileNotFoundError(f"Edge cookie database not found: {cookies_db}")

    master_key = get_edge_master_key()

    # Copy DB to temp location using Windows API with FILE_SHARE_READ|WRITE|DELETE
    # (bypasses the exclusive lock Edge holds via normal Python open)
    tmp_db = os.path.join(tempfile.gettempdir(), 'edge_cookies_tmp.db')
    import ctypes
    import ctypes.wintypes as wt

    GENERIC_READ = 0x80000000
    FILE_SHARE_READ  = 0x00000001
    FILE_SHARE_WRITE = 0x00000002
    FILE_SHARE_DELETE = 0x00000004
    OPEN_EXISTING = 3
    FILE_FLAG_SEQUENTIAL_SCAN = 0x08000000

    kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
    h = kernel32.CreateFileW(
        cookies_db,
        GENERIC_READ,
        FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
        None, OPEN_EXISTING,
        FILE_FLAG_SEQUENTIAL_SCAN, None
    )
    if h == wt.HANDLE(-1).value:
        raise OSError(f"Cannot open cookie DB: error {ctypes.get_last_error()}")
    try:
        file_size = wt.LARGE_INTEGER(0)
        kernel32.GetFileSizeEx(h, ctypes.byref(file_size))
        size = file_size.value
        buf = (ctypes.c_char * size)()
        read = wt.DWORD(0)
        kernel32.ReadFile(h, buf, size, ctypes.byref(read), None)
        raw_data = bytes(buf[:read.value])
    finally:
        kernel32.CloseHandle(h)

    with open(tmp_db, 'wb') as f:
        f.write(raw_data)

    result = {}
    try:
        conn = sqlite3.connect(tmp_db)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT name, encrypted_value FROM cookies WHERE host_key LIKE '%x.com%' AND name IN ('auth_token','ct0')"
        )
        for name, enc_val in cursor.fetchall():
            value = decrypt_cookie_value(enc_val, master_key)
            if value:
                result[name] = value
        conn.close()
    finally:
        os.remove(tmp_db)

    return result


def scrape_grok_share(url):
    """Fetch Grok share URL with authentication and extract conversation content."""
    print(f"Extracting Edge cookies for x.com...")
    cookies = get_xcom_cookies()

    if not cookies:
        print("ERROR: No x.com auth cookies found. Are you logged into X in Edge?")
        sys.exit(1)

    found = [k for k in ['auth_token', 'ct0'] if k in cookies]
    print(f"  Found cookies: {found}")

    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36 Edg/120.0.0.0',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Referer': 'https://x.com/',
    })
    for name, value in cookies.items():
        session.cookies.set(name, value, domain='.x.com')

    print(f"Fetching: {url}")
    resp = session.get(url, timeout=30)
    print(f"  HTTP {resp.status_code}, size: {len(resp.content)} bytes")

    html = resp.text

    # --- Extract OG meta content (server-side rendered content block) ---
    og_title = ''
    og_desc = ''
    m = re.search(r'<meta[^>]*property=["\']og:title["\'][^>]*content=["\']([^"\']+)["\']', html)
    if m:
        og_title = m.group(1)
    m = re.search(r'<meta[^>]*property=["\']og:description["\'][^>]*content=["\']([^"\']+)["\']', html, re.DOTALL)
    if not m:
        # Try flipped attribute order
        m = re.search(r'<meta[^>]*content=["\']([^"\']{50,})["\'][^>]*property=["\']og:description["\']', html, re.DOTALL)
    if m:
        og_desc = m.group(1)

    # --- Broader content extraction: look for all large meta content blocks ---
    all_meta = re.findall(r'<meta content="([\s\S]{100,?}?)"[^>]*data-rh', html)

    content_parts = []
    if og_title:
        content_parts.append(f"[OG TITLE]: {og_title}")
    if og_desc:
        content_parts.append(f"[OG DESCRIPTION]:\n{og_desc}")
    for block in all_meta:
        if block not in [og_title, og_desc]:
            content_parts.append(f"[META BLOCK]:\n{block}")

    if not content_parts:
        # Fallback: dump all large text content from HTML
        # Remove script/style tags
        clean = re.sub(r'<script[^>]*>[\s\S]*?</script>', '', html)
        clean = re.sub(r'<style[^>]*>[\s\S]*?</style>', '', clean)
        # Strip tags
        clean = re.sub(r'<[^>]+>', ' ', clean)
        clean = re.sub(r'\s+', ' ', clean).strip()
        content_parts.append(f"[PAGE TEXT]:\n{clean[:50000]}")

    full_content = '\n\n'.join(content_parts)

    # Extract ID for filenames
    url_id_match = re.search(r'/share/([a-f0-9A-F]+)', url)
    url_id = url_id_match.group(1) if url_id_match else 'unknown'

    content_file = f'grok_share_{url_id}_content.txt'
    html_file = f'grok_share_{url_id}_source.html'

    with open(content_file, 'w', encoding='utf-8') as f:
        f.write(full_content)
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"\n{'='*80}")
    print("EXTRACTED CONTENT:")
    print('='*80)
    print(full_content[:3000])
    if len(full_content) > 3000:
        print(f"... [{len(full_content) - 3000} more chars] ...")
    print('='*80)
    print(f"\n✓ Content saved to: {content_file} ({len(full_content)} chars)")
    print(f"✓ HTML source saved to: {html_file} ({len(html)} chars)")

    return full_content


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python scrape_grok_auth.py <URL>")
        sys.exit(1)
    scrape_grok_share(sys.argv[1])
