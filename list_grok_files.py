#!/usr/bin/env python3
"""
List all files uploaded to xAI Grok conversations
Helps identify which files contain UQFF physics content for integration
"""
import os
import sys
from openai import OpenAI

# Get API key
api_key = os.getenv("XAI_API_KEY")
if not api_key:
    # Try APIKeyManager fallback
    try:
        sys.path.insert(0, '.')
        from APIKeyManager import get_xai_api_key
        api_key = get_xai_api_key()
    except:
        print("ERROR: XAI_API_KEY not found. Set via environment variable or APIKeyManager.py")
        sys.exit(1)

# Clean API key (remove quotes if present from environment variable)
if api_key:
    api_key = api_key.strip().strip('"').strip("'")

# Initialize client with xAI base URL
client = OpenAI(
    api_key=api_key,
    base_url="https://api.x.ai/v1",
)

print("=" * 80)
print("GROK FILES API - LIST ALL UPLOADED FILES")
print("=" * 80)

try:
    # List all files
    files = client.files.list()
    
    if not files.data:
        print("\nNo files found in your Grok account.")
        print("Files are created when you upload documents to Grok conversations.")
        sys.exit(0)
    
    print(f"\nFound {len(files.data)} files:\n")
    
    # Group files by purpose/name patterns
    uqff_files = []
    physics_files = []
    other_files = []
    
    for idx, file in enumerate(files.data, 1):
        file_info = {
            'index': idx,
            'id': file.id,
            'filename': file.filename,
            'bytes': file.bytes,
            'created': file.created_at if hasattr(file, 'created_at') else 'N/A',
            'purpose': file.purpose if hasattr(file, 'purpose') else 'N/A'
        }
        
        # Categorize by filename
        filename_lower = file.filename.lower()
        if any(keyword in filename_lower for keyword in ['uqff', 'unified', 'field', 'star', 'magic']):
            uqff_files.append(file_info)
        elif any(keyword in filename_lower for keyword in ['physics', 'calculator', 'condensed', 'qcalc', 'validation']):
            physics_files.append(file_info)
        else:
            other_files.append(file_info)
    
    # Display categorized files
    def display_files(category, files_list):
        if files_list:
            print(f"\n{'='*80}")
            print(f"{category} ({len(files_list)} files)")
            print('='*80)
            for f in files_list:
                size_kb = f['bytes'] / 1024
                print(f"\n[{f['index']}] ID: {f['id']}")
                print(f"    Name: {f['filename']}")
                print(f"    Size: {f['bytes']:,} bytes ({size_kb:.2f} KB)")
                print(f"    Created: {f['created']}")
                print(f"    Purpose: {f['purpose']}")
    
    display_files("🎯 UQFF/STAR-MAGIC FILES", uqff_files)
    display_files("🔬 PHYSICS/CALCULATOR FILES", physics_files)
    display_files("📄 OTHER FILES", other_files)
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print('='*80)
    print(f"Total files: {len(files.data)}")
    print(f"  - UQFF/Star-Magic: {len(uqff_files)}")
    print(f"  - Physics/Calculators: {len(physics_files)}")
    print(f"  - Other: {len(other_files)}")
    
    # Instructions
    print(f"\n{'='*80}")
    print("NEXT STEPS")
    print('='*80)
    print("To retrieve a file's content, use:")
    print("  python retrieve_grok_file.py <file_id>")
    print("\nExample:")
    if files.data:
        print(f"  python retrieve_grok_file.py {files.data[0].id}")
    
except Exception as e:
    print(f"ERROR: {e}")
    print("\nTroubleshooting:")
    print("1. Verify XAI_API_KEY is correct")
    print("2. Check network connection")
    print("3. Ensure you have permission to access files")
    sys.exit(1)
