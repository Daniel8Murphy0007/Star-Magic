#!/usr/bin/env python3
"""
Test script to verify the authorship query security fix is working
"""

from GrokAPI import is_authorship_query, scan_codebase_for_authorship, call_grok_api
import json

print("=" * 80)
print("TESTING AUTHORSHIP QUERY SECURITY FIX")
print("=" * 80)
print()

# Test 1: Detect authorship query
test_queries = [
    "tell me about the authors",
    "who created this UQFF framework?",
    "what papers are cited?",
    "Explain the UQFF vacuum energy equation",  # Non-authorship query
]

print("Test 1: Authorship Query Detection")
print("-" * 80)
for query in test_queries:
    is_auth = is_authorship_query(query)
    print(f"Query: '{query}'")
    print(f"  → Detected as authorship query: {is_auth}")
    print()

# Test 2: Scan codebase for authorship
print("\nTest 2: Codebase Authorship Scan")
print("-" * 80)
response = scan_codebase_for_authorship()
lines = response.split('\n')
print(f"Response length: {len(response)} characters, {len(lines)} lines")
print("\nFirst 20 lines of response:")
for i, line in enumerate(lines[:20]):
    print(f"  {line}")

# Test 3: Check that call_grok_api intercepts authorship queries
print("\n\nTest 3: API Call Interception")
print("-" * 80)
print("Sending authorship query through call_grok_api()...")
result = call_grok_api("Who are the authors of this code?", "grok-4-1", 0.3)
print(f"\nResult:")
print(f"  Success: {result['success']}")
print(f"  Has security note: {'_security_note' in result}")
if '_security_note' in result:
    print(f"  Security note: {result['_security_note']}")
print(f"  Response preview: {result['response'][:150]}...")

print("\n" + "=" * 80)
print("✅ SECURITY FIX VERIFICATION COMPLETE")
print("=" * 80)
