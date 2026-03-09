#!/usr/bin/env python3
"""
Retrieve and download a specific file from xAI Grok by file ID
Usage: python retrieve_grok_file.py <file_id>
"""
import os
import sys
import json
from datetime import datetime
from openai import OpenAI

if len(sys.argv) < 2:
    print("Usage: python retrieve_grok_file.py <file_id>")
    print("\nTo find file IDs, run: python list_grok_files.py")
    sys.exit(1)

file_id = sys.argv[1]

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

# Initialize client
client = OpenAI(
    api_key=api_key,
    base_url="https://api.x.ai/v1",
)

print("=" * 80)
print(f"RETRIEVING FILE: {file_id}")
print("=" * 80)

try:
    # Get file metadata
    print("\n[1/3] Fetching metadata...")
    file = client.files.retrieve(file_id)
    
    print(f"\nFilename: {file.filename}")
    print(f"Size: {file.bytes:,} bytes ({file.bytes/1024:.2f} KB)")
    print(f"Purpose: {file.purpose if hasattr(file, 'purpose') else 'N/A'}")
    
    # Get file content
    print("\n[2/3] Downloading content...")
    content_response = client.files.content(file_id)
    content = content_response.content  # As bytes
    
    print(f"Downloaded: {len(content):,} bytes")
    
    # Determine output filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_name = file.filename.rsplit('.', 1)[0] if '.' in file.filename else file.filename
    extension = file.filename.rsplit('.', 1)[1] if '.' in file.filename else 'txt'
    output_file = f"grok_file_{base_name}_{timestamp}.{extension}"
    
    # Save content
    print(f"\n[3/3] Saving to: {output_file}")
    with open(output_file, "wb") as f:
        f.write(content)
    
    # Try to decode as text for preview
    try:
        text_content = content.decode("utf-8")
        print(f"\n{'='*80}")
        print("CONTENT PREVIEW (first 2000 chars)")
        print('='*80)
        print(text_content[:2000])
        if len(text_content) > 2000:
            print(f"\n...{len(text_content)-2000:,} more characters...")
        
        # Analyze content for UQFF keywords
        print(f"\n{'='*80}")
        print("CONTENT ANALYSIS")
        print('='*80)
        
        keywords = {
            'UQFF': ['uqff', 'unified quantum field', 'unified field'],
            'Equations': ['f_u', 'ug1', 'ug2', 'ug3', 'ug4', 'f_ubi', 'f_u_bi_i'],
            'Constants': ['kappa', 'ssq', 'h_scm', 'u_ua', 'beta_i'],
            'Systems': ['sgr1745', 'm87', 'sagittarius', 'magnetar', 'black hole'],
            'Calculators': ['calculator', 'class', 'def compute', 'condensedphysics'],
            'Validation': ['validation', 'test', 'pass', 'fail', 'accuracy']
        }
        
        text_lower = text_content.lower()
        found_keywords = {}
        for category, terms in keywords.items():
            count = sum(text_lower.count(term) for term in terms)
            if count > 0:
                found_keywords[category] = count
        
        if found_keywords:
            print("Found UQFF-related keywords:")
            for category, count in sorted(found_keywords.items(), key=lambda x: x[1], reverse=True):
                print(f"  - {category}: {count} occurrences")
        else:
            print("⚠️  No UQFF-specific keywords found")
            print("   This file may not contain Star-Magic physics content")
        
        # Check for Python classes
        if 'class ' in text_content:
            import re
            classes = re.findall(r'class\s+(\w+)\s*[:\(]', text_content)
            if classes:
                print(f"\nPython classes found ({len(classes)}):")
                for cls in classes[:20]:  # First 20
                    print(f"  - {cls}")
                if len(classes) > 20:
                    print(f"  ... and {len(classes)-20} more")
        
    except UnicodeDecodeError:
        print("\n⚠️  Content is binary (PDF, image, etc.) - cannot preview as text")
    
    print(f"\n{'='*80}")
    print("SUCCESS")
    print('='*80)
    print(f"✅ File saved: {output_file}")
    print(f"✅ Size: {len(content):,} bytes")
    
    # Next steps
    print(f"\n{'='*80}")
    print("NEXT STEPS FOR INTEGRATION")
    print('='*80)
    print("1. Review the file content above")
    print("2. If it contains UQFF physics:")
    print("   - Extract calculator classes")
    print("   - Check for duplicates: grep search vs CondensedPhysics.py (81K lines)")
    print("   - Add unique calculators to CondensedPhysics.py")
    print("3. If it's validation data:")
    print("   - Check file format (CSV, JSON, etc.)")
    print("   - Add to validation pipeline")
    
except Exception as e:
    print(f"\nERROR: {e}")
    print("\nPossible causes:")
    print("- Invalid file ID (run list_grok_files.py to see available files)")
    print("- File was deleted")
    print("- Network error")
    sys.exit(1)
