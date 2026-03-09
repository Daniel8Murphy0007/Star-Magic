#!/usr/bin/env python3
"""
Download ALL files from xAI Grok account for offline analysis
Saves files with timestamps and creates an index
"""
import os
import sys
import json
from datetime import datetime
from openai import OpenAI

# Get API key
api_key = os.getenv("XAI_API_KEY")
if not api_key:
    try:
        sys.path.insert(0, '.')
        from APIKeyManager import get_xai_api_key
        api_key = get_xai_api_key()
    except:
        print("ERROR: XAI_API_KEY not found")
        sys.exit(1)

# Clean API key (remove quotes if present from environment variable)
if api_key:
    api_key = api_key.strip().strip('"').strip("'")

# Create output directory
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_dir = f"grok_files_backup_{timestamp}"
os.makedirs(output_dir, exist_ok=True)

print("=" * 80)
print(f"DOWNLOADING ALL GROK FILES TO: {output_dir}/")
print("=" * 80)

# Initialize client
client = OpenAI(
    api_key=api_key,
    base_url="https://api.x.ai/v1",
)

try:
    # List all files
    print("\n[1/3] Fetching file list...")
    files = client.files.list()
    
    if not files.data:
        print("No files found.")
        sys.exit(0)
    
    print(f"Found {len(files.data)} files\n")
    
    # Download each file
    print("[2/3] Downloading files...")
    manifest = []
    
    for idx, file in enumerate(files.data, 1):
        print(f"\n[{idx}/{len(files.data)}] {file.filename} ({file.bytes:,} bytes)")
        
        try:
            # Get content
            content_response = client.files.content(file.id)
            content = content_response.content
            
            # Save with safe filename
            safe_name = "".join(c if c.isalnum() or c in "._- " else "_" for c in file.filename)
            output_path = os.path.join(output_dir, f"{idx:03d}_{safe_name}")
            
            with open(output_path, "wb") as f:
                f.write(content)
            
            print(f"  ✅ Saved: {output_path}")
            
            # Add to manifest
            file_info = {
                'index': idx,
                'id': file.id,
                'filename': file.filename,
                'saved_as': output_path,
                'size_bytes': file.bytes,
                'created_at': str(file.created_at) if hasattr(file, 'created_at') else 'N/A',
                'purpose': file.purpose if hasattr(file, 'purpose') else 'N/A'
            }
            
            # Try to extract summary
            try:
                text = content.decode('utf-8')
                lines = text.split('\n')
                file_info['preview'] = ' '.join(lines[:5])[:200]
                file_info['line_count'] = len(lines)
                
                # Check for UQFF
                text_lower = text.lower()
                file_info['contains_uqff'] = any(term in text_lower for term in ['uqff', 'unified field', 'f_u', 'ug1'])
                file_info['contains_calculators'] = 'class ' in text and 'Calculator' in text
                
            except:
                file_info['preview'] = '[Binary file]'
                file_info['contains_uqff'] = False
                file_info['contains_calculators'] = False
            
            manifest.append(file_info)
            
        except Exception as e:
            print(f"  ❌ ERROR: {e}")
    
    # Save manifest
    print(f"\n[3/3] Creating manifest...")
    manifest_path = os.path.join(output_dir, "MANIFEST.json")
    with open(manifest_path, "w", encoding='utf-8') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"✅ Manifest saved: {manifest_path}")
    
    # Create readable index
    index_path = os.path.join(output_dir, "INDEX.txt")
    with open(index_path, "w", encoding='utf-8') as f:
        f.write("GROK FILES BACKUP INDEX\n")
        f.write("=" * 80 + "\n")
        f.write(f"Downloaded: {datetime.now()}\n")
        f.write(f"Total files: {len(manifest)}\n\n")
        
        # Group by UQFF relevance
        uqff_files = [m for m in manifest if m.get('contains_uqff')]
        calc_files = [m for m in manifest if m.get('contains_calculators')]
        
        if uqff_files:
            f.write("\n🎯 UQFF-RELATED FILES:\n")
            f.write("-" * 80 + "\n")
            for m in uqff_files:
                f.write(f"\n[{m['index']}] {m['filename']}\n")
                f.write(f"    File: {m['saved_as']}\n")
                f.write(f"    Size: {m['size_bytes']:,} bytes\n")
                f.write(f"    Preview: {m.get('preview', 'N/A')}\n")
        
        if calc_files:
            f.write("\n\n🔬 CALCULATOR FILES:\n")
            f.write("-" * 80 + "\n")
            for m in calc_files:
                f.write(f"\n[{m['index']}] {m['filename']}\n")
                f.write(f"    File: {m['saved_as']}\n")
                f.write(f"    Size: {m['size_bytes']:,} bytes\n")
        
        f.write("\n\nALL FILES:\n")
        f.write("-" * 80 + "\n")
        for m in manifest:
            f.write(f"[{m['index']:03d}] {m['filename']:50s} {m['size_bytes']:>10,} bytes\n")
    
    print(f"✅ Index saved: {index_path}")
    
    # Summary
    print(f"\n{'='*80}")
    print("DOWNLOAD COMPLETE")
    print('='*80)
    print(f"📁 Directory: {output_dir}/")
    print(f"📄 Files downloaded: {len(manifest)}")
    print(f"🎯 UQFF-related: {len(uqff_files)}")
    print(f"🔬 Calculator files: {len(calc_files)}")
    print(f"\nNext steps:")
    print(f"1. Review INDEX.txt for UQFF-related files")
    print(f"2. Open downloaded files for physics extraction")
    print(f"3. Run grep searches to check for duplicates")
    
except Exception as e:
    print(f"\nERROR: {e}")
    sys.exit(1)
