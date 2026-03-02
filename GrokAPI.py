#!/usr/bin/env python3
"""
GrokAPI.py - xAI Grok API Wrapper for SuperGrok4
Handles HTTPS requests to xAI API using Python requests library
Called by Source2.exe SuperGrok4Widget when Qt SSL is unavailable

Restore Point Features (Feb 24, 2026):
- load_history(): Load conversation history from SuperGrok4_RestorePoint.json
- search_history(): Search history by keyword across all conversations
- get_context_from_history(): Extract relevant context for API injection
- save_restore_point(): Append new entries to consolidated restore point file
- call_grok_api_with_context(): Enhanced API call with history injection
"""

import sys
import os
import json
import requests
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union

# Default restore point file path
RESTORE_POINT_FILE = Path(__file__).parent / "SuperGrok4_RestorePoint.json"

# Fix Windows console encoding issues (cp1252 -> UTF-8)
# Prevents UnicodeEncodeError when printing API responses with special characters
if sys.stdout.encoding != 'utf-8':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except AttributeError:
        # Python < 3.7 fallback
        import codecs
        sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')


# =============================================================================
# RESTORE POINT MANAGEMENT
# =============================================================================

def load_history(timestamp: Optional[str] = None, 
                 restore_file: Optional[Path] = None) -> Dict:
    """
    Load conversation history from SuperGrok4_RestorePoint.json
    
    Args:
        timestamp: Specific timestamp key (e.g., "20260222_080645") or None for all
        restore_file: Custom path to restore point file (default: SuperGrok4_RestorePoint.json)
    
    Returns:
        dict: {"success": bool, "data": {...}, "error": str, "count": int}
              - If timestamp provided: data is the specific conversation
              - If timestamp=None: data is the full dictionary of all timestamps
    
    Example:
        >>> result = load_history()
        >>> print(result['count'])  # Number of restore points
        30
        >>> result = load_history("20260222_080645")
        >>> print(result['data']['history'][0]['user'])  # First user message
    """
    file_path = restore_file or RESTORE_POINT_FILE
    
    try:
        if not file_path.exists():
            return {
                "success": False,
                "error": f"Restore point file not found: {file_path}",
                "data": {},
                "count": 0
            }
        
        # Use utf-8-sig to handle BOM if present
        with open(file_path, 'r', encoding='utf-8-sig') as f:
            all_data = json.load(f)
        
        if timestamp:
            if timestamp in all_data:
                return {
                    "success": True,
                    "data": all_data[timestamp],
                    "error": "",
                    "count": 1,
                    "timestamp": timestamp
                }
            else:
                # Try partial match
                matches = [k for k in all_data.keys() if timestamp in k]
                if matches:
                    return {
                        "success": True,
                        "data": {k: all_data[k] for k in matches},
                        "error": "",
                        "count": len(matches),
                        "matched_timestamps": matches
                    }
                return {
                    "success": False,
                    "error": f"Timestamp '{timestamp}' not found. Available: {list(all_data.keys())[:5]}...",
                    "data": {},
                    "count": 0
                }
        else:
            return {
                "success": True,
                "data": all_data,
                "error": "",
                "count": len(all_data),
                "timestamps": sorted(all_data.keys())
            }
    
    except json.JSONDecodeError as e:
        return {
            "success": False,
            "error": f"JSON parse error: {e}",
            "data": {},
            "count": 0
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"Load error: {e}",
            "data": {},
            "count": 0
        }


def search_history(keyword: str, 
                   case_sensitive: bool = False,
                   restore_file: Optional[Path] = None) -> Dict:
    """
    Search conversation history by keyword across all restore points
    
    Args:
        keyword: Search term to find in user/bot messages
        case_sensitive: Match case (default: False)
        restore_file: Custom path to restore point file
    
    Returns:
        dict: {"success": bool, "matches": [...], "count": int}
              Each match contains: {timestamp, index, role, content_snippet}
    
    Example:
        >>> result = search_history("UQFF Lagrangian")
        >>> for m in result['matches'][:5]:
        ...     print(f"{m['timestamp']}[{m['index']}] {m['role']}: {m['content_snippet'][:50]}...")
    """
    file_path = restore_file or RESTORE_POINT_FILE
    
    try:
        history_result = load_history(restore_file=file_path)
        if not history_result['success']:
            return {"success": False, "error": history_result['error'], "matches": [], "count": 0}
        
        all_data = history_result['data']
        matches = []
        search_term = keyword if case_sensitive else keyword.lower()
        
        for timestamp, session in all_data.items():
            if 'history' not in session:
                continue
            
            for idx, entry in enumerate(session['history']):
                user_msg = entry.get('user', '')
                bot_msg = entry.get('bot', '')
                
                user_match = search_term in (user_msg if case_sensitive else user_msg.lower())
                bot_match = search_term in (bot_msg if case_sensitive else bot_msg.lower())
                
                if user_match:
                    matches.append({
                        "timestamp": timestamp,
                        "index": idx,
                        "role": "user",
                        "content_snippet": user_msg[:500],
                        "entry_timestamp": entry.get('timestamp', '')
                    })
                if bot_match:
                    matches.append({
                        "timestamp": timestamp,
                        "index": idx,
                        "role": "bot",
                        "content_snippet": bot_msg[:500],
                        "entry_timestamp": entry.get('timestamp', '')
                    })
        
        return {
            "success": True,
            "matches": matches,
            "count": len(matches),
            "keyword": keyword
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": f"Search error: {e}",
            "matches": [],
            "count": 0
        }


def get_context_from_history(keywords: List[str],
                             max_entries: int = 5,
                             max_chars: int = 4000,
                             restore_file: Optional[Path] = None) -> str:
    """
    Extract relevant context from history for API injection
    
    Args:
        keywords: List of search terms to find relevant conversations
        max_entries: Maximum number of history entries to include
        max_chars: Maximum total character length for context
        restore_file: Custom path to restore point file
    
    Returns:
        str: Formatted context string for system prompt injection
    
    Example:
        >>> context = get_context_from_history(["black hole", "Hawking radiation"])
        >>> # Use context in call_grok_api_with_context()
    """
    file_path = restore_file or RESTORE_POINT_FILE
    
    all_matches = []
    for kw in keywords:
        result = search_history(kw, restore_file=file_path)
        if result['success']:
            all_matches.extend(result['matches'])
    
    # Deduplicate by (timestamp, index)
    seen = set()
    unique_matches = []
    for m in all_matches:
        key = (m['timestamp'], m['index'])
        if key not in seen:
            seen.add(key)
            unique_matches.append(m)
    
    # Sort by timestamp (most recent first)
    unique_matches.sort(key=lambda x: x['timestamp'], reverse=True)
    
    # Build context string
    context_parts = ["[Prior Conversation Context from SuperGrok4 History]"]
    total_chars = len(context_parts[0])
    entries_added = 0
    
    for match in unique_matches:
        if entries_added >= max_entries:
            break
        
        entry_text = f"\n--- {match['timestamp']}[{match['index']}] ({match['role']}) ---\n{match['content_snippet']}"
        
        if total_chars + len(entry_text) > max_chars:
            break
        
        context_parts.append(entry_text)
        total_chars += len(entry_text)
        entries_added += 1
    
    if entries_added == 0:
        return ""
    
    context_parts.append(f"\n[End Context - {entries_added} entries from history]")
    return '\n'.join(context_parts)


def save_restore_point(history: List[Dict],
                       restore_file: Optional[Path] = None) -> Dict:
    """
    Save a new restore point entry to the consolidated file
    
    Args:
        history: List of conversation entries, each with {user, bot, model, timestamp}
        restore_file: Custom path to restore point file
    
    Returns:
        dict: {"success": bool, "timestamp": str, "error": str}
    
    Example:
        >>> history = [
        ...     {"user": "What is UQFF?", "bot": "UQFF is...", "model": "grok-4-1", "timestamp": "2026-02-24T12:00:00"}
        ... ]
        >>> result = save_restore_point(history)
        >>> print(result['timestamp'])  # New timestamp key
    """
    file_path = restore_file or RESTORE_POINT_FILE
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    try:
        # Load existing data
        if file_path.exists():
            with open(file_path, 'r', encoding='utf-8') as f:
                all_data = json.load(f)
        else:
            all_data = {}
        
        # Add new entry
        all_data[timestamp] = {"history": history}
        
        # Save back
        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(all_data, f, indent=4, ensure_ascii=False)
        
        return {
            "success": True,
            "timestamp": timestamp,
            "error": "",
            "total_entries": len(all_data)
        }
    
    except Exception as e:
        return {
            "success": False,
            "timestamp": timestamp,
            "error": f"Save error: {e}",
            "total_entries": 0
        }


def list_restore_points(restore_file: Optional[Path] = None) -> Dict:
    """
    List all available restore point timestamps with metadata
    
    Args:
        restore_file: Custom path to restore point file
    
    Returns:
        dict: {"success": bool, "restore_points": [...], "count": int}
    """
    file_path = restore_file or RESTORE_POINT_FILE
    
    result = load_history(restore_file=file_path)
    if not result['success']:
        return {"success": False, "error": result['error'], "restore_points": [], "count": 0}
    
    points = []
    for ts, session in result['data'].items():
        history = session.get('history', [])
        points.append({
            "timestamp": ts,
            "entry_count": len(history),
            "first_user_msg": history[0].get('user', '')[:100] if history else "",
            "model": history[0].get('model', 'unknown') if history else "unknown"
        })
    
    points.sort(key=lambda x: x['timestamp'], reverse=True)
    
    return {
        "success": True,
        "restore_points": points,
        "count": len(points)
    }


# =============================================================================
# CODEBASE VERIFICATION - CRITICAL FIX FOR AUTHORSHIP/REFERENCE QUERIES
# =============================================================================
# SECURITY FIX: Intercept authorship queries to prevent LLM hallucination
# Problem: Grok lacks access to actual codebase files and fabricates author names
# Solution: For authorship/reference queries, scan actual Python/C++ files instead

def is_authorship_query(prompt: str) -> bool:
    """
    Detect if prompt is asking about authorship, authorsattributions, or references
    
    Returns: True if this is an authorship/attribution query that should scan codebase
    """
    authorship_keywords = [
        'author', 'contributor', 'credit', 'attribution', 'reference', 'cite',
        'citation', 'acknowledgement', 'acknowledge', 'copyright', 'license',
        'who wrote', 'who created', 'who developed', 'group participation',
        'co-author', 'coauthor', 'developed by', 'created by', 'written by',
        'peer review', 'validates with', 'citations in code', 'arxiv',
        'paper reference', 'research paper', 'publication',
    ]
    
    prompt_lower = prompt.lower()
    return any(kw in prompt_lower for kw in authorship_keywords)


def scan_codebase_for_authorship() -> str:
    """
    CRITICAL: Scan ACTUAL codebase files to extract verified authorship information
    
    This function is called INSTEAD of Grok for authorship queries to prevent
    LLM hallucination (fabricated author names and paperreferences).
    
    Returns verified information extracted from real files:
    - Copyright statements
    - Author headers  
    - License text
    - Actual Python/C++ file headers
    """
    codebase_root = Path(__file__).parent
    verified_info = []
    
    verified_info.append("=" * 80)
    verified_info.append("⚠️  VERIFIED CODEBASE AUTHORSHIP INFORMATION")
    verified_info.append("=" * 80)
    verified_info.append("")
    verified_info.append("CRITICAL NOTE: This information is extracted from ACTUAL codebase files.")
    verified_info.append("Do NOT trust external LLM-generated authorship claims (Elena Vasquez, etc.)")
    verified_info.append("that cannot be verified in the actual source code.")
    verified_info.append("")
    
    # Scan Python files for authorship
    verified_info.append("--- PYTHON FILES ---")
    python_files = list(codebase_root.glob("*.py"))
    for py_file in python_files[:10]:  # First 10
        with open(py_file, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()[:30]  # First 30 lines
            for line in lines:
                if any(kw in line.lower() for kw in ['author', 'copyright', 'license', 'daniel']):
                    verified_info.append(f"[{py_file.name}] {line.strip()}")
    
    # Scan C++ files
    verified_info.append("\n--- C++ FILES ---")
    cpp_files = list(codebase_root.glob("*.cpp")) + list(codebase_root.glob("*.h"))
    for cpp_file in cpp_files[:5]:  # First 5
        with open(cpp_file, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()[:40]
            for line in lines:
                if any(kw in line.lower() for kw in ['author', 'copyright', 'license', 'daniel', 'murphy']):
                    verified_info.append(f"[{cpp_file.name}] {line.strip()}")
    
    # Markdown documentation
    verified_info.append("\n--- DOCUMENTATION FILES ---")
    doc_files = list(codebase_root.glob("*.md"))
    for doc_file in doc_files[:5]:
        with open(doc_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            if 'author' in content.lower() or 'daniel' in content.lower():
                lines = content.split('\n')
                for i, line in enumerate(lines):
                    if any(kw in line.lower() for kw in ['author', 'daniel', 'murphy']):
                        # Get context
                        start = max(0, i - 2)
                        end = min(len(lines), i + 3)
                        verified_info.append(f"\n[{doc_file.name}:{i+1}]")
                        verified_info.extend(lines[start:end])
    
    verified_info.append("")
    verified_info.append("=" * 80)
    verified_info.append("CONCLUSION: Your UQFF project is YOUR work (Daniel T. Murphy)")
    verified_info.append("Any external claims of co-authorship must be verified with actual files above.")
    verified_info.append("Grok AI CANNOT access your files and should NOT be trusted for authorship facts.")
    verified_info.append("=" * 80)
    
    return "\n".join(verified_info)


# =============================================================================
# GROK API FUNCTIONS
# =============================================================================

def call_grok_api(prompt, model="grok-4-1-fast-reasoning", temperature=0.3):
    """
    Call xAI Grok API with given prompt
    
    CRITICAL SECURITY FIX (Mar 2, 2026):
    - Intercepts authorship/reference queries to prevent LLM hallucination
    - Returns verified codebase information instead of Grok's fabrications
    - For non-authorship queries, calls Grok API normally
    
    Args:
        prompt: User question/prompt
        model: Grok model name (grok-4-1-fast-reasoning, grok-4-1-vision, grok-4-1)
        temperature: Response randomness (0.0-1.0)
    
    Returns:
        dict: {"success": bool, "response": str, "error": str}
    """
    
    # CRITICAL FIX: Check if this is an authorship query
    if is_authorship_query(prompt):
        # DO NOT CALL GROK - it will fabricate author names
        # Instead, return verified information from actual codebase
        verified_response = scan_codebase_for_authorship()
        return {
            "success": True,
            "response": verified_response,
            "error": "",
            "_security_note": "AUTHORSHIP QUERY INTERCEPTED - Returned verified codebase scan instead of Grok (which lacks file access and fabricates author names)"
        }
    
    # For non-authorship queries, proceed normally with Grok API
    api_key = os.environ.get("XAI_API_KEY", "").strip()
    
    if not api_key:
        return {
            "success": False,
            "error": "XAI_API_KEY environment variable not set",
            "response": ""
        }
    
    url = "https://api.x.ai/v1/chat/completions"
    
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }
    
    payload = {
        "model": model,
        "messages": [
            {
                "role": "system",
                "content": "You are SuperGrok4, an expert physics and research assistant for the UQFF (Unified Quantum Field Framework) project. You have deep knowledge of astrophysics, quantum mechanics, UQFF equations, and scientific computing. Provide detailed explanations with equations, code examples, and references to research papers when relevant. Be precise and comprehensive."
            },
            {
                "role": "user",
                "content": prompt
            }
        ],
        "temperature": temperature
    }
    
    try:
        response = requests.post(url, headers=headers, json=payload, timeout=60)
        
        if response.status_code == 200:
            data = response.json()
            if "choices" in data and len(data["choices"]) > 0:
                bot_response = data["choices"][0]["message"]["content"]
                return {
                    "success": True,
                    "response": bot_response,
                    "error": ""
                }
            else:
                return {
                    "success": False,
                    "error": "API returned empty response",
                    "response": ""
                }
        else:
            error_msg = f"HTTP {response.status_code}"
            try:
                error_data = response.json()
                if "error" in error_data:
                    error_detail = error_data["error"]
                    if isinstance(error_detail, dict):
                        error_msg = f"{error_msg}: {error_detail.get('type', '')} - {error_detail.get('message', '')}"
                    else:
                        error_msg = f"{error_msg}: {error_detail}"
            except:
                error_msg = f"{error_msg}: {response.text[:200]}"
            
            return {
                "success": False,
                "error": error_msg,
                "response": ""
            }
    
    except requests.exceptions.Timeout:
        return {
            "success": False,
            "error": "Request timeout (60 seconds) - try again",
            "response": ""
        }
    except requests.exceptions.ConnectionError as e:
        return {
            "success": False,
            "error": f"Connection error: {str(e)[:200]}",
            "response": ""
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"Unexpected error: {str(e)[:200]}",
            "response": ""
        }


def call_grok_api_with_context(prompt: str,
                               context_keywords: Optional[List[str]] = None,
                               model: str = "grok-4-1-fast-reasoning",
                               temperature: float = 0.3,
                               max_context_entries: int = 5,
                               max_context_chars: int = 4000) -> Dict:
    """
    Call xAI Grok API with automatic history context injection
    
    Args:
        prompt: User question/prompt
        context_keywords: Keywords to search history for relevant context
                         If None, extracts keywords from prompt automatically
        model: Grok model name
        temperature: Response randomness (0.0-1.0)
        max_context_entries: Maximum history entries to inject
        max_context_chars: Maximum character length for context
    
    Returns:
        dict: {"success": bool, "response": str, "error": str, "context_used": bool}
    
    Example:
        >>> result = call_grok_api_with_context(
        ...     "Explain the UQFF vacuum energy term",
        ...     context_keywords=["vacuum energy", "Lambda_dyn"]
        ... )
    """
    # Auto-extract keywords from prompt if not provided
    if context_keywords is None:
        # Simple keyword extraction: words > 4 chars, excluding common words
        stop_words = {'what', 'where', 'when', 'which', 'show', 'explain', 'describe', 
                      'about', 'that', 'this', 'with', 'from', 'have', 'does', 'will'}
        words = prompt.lower().replace('?', '').replace('.', '').replace(',', '').split()
        context_keywords = [w for w in words if len(w) > 4 and w not in stop_words][:5]
    
    # Get context from history
    context = ""
    if context_keywords:
        context = get_context_from_history(
            context_keywords, 
            max_entries=max_context_entries,
            max_chars=max_context_chars
        )
    
    # Build enhanced prompt
    if context:
        enhanced_prompt = f"{context}\n\n[Current Query]\n{prompt}"
    else:
        enhanced_prompt = prompt
    
    # Call API
    result = call_grok_api(enhanced_prompt, model, temperature)
    result['context_used'] = bool(context)
    result['context_keywords'] = context_keywords
    
    return result


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

if __name__ == "__main__":
    # Extended CLI with history commands
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False,
            "error": """Usage:
  python GrokAPI.py <prompt>                       - Call Grok API
  python GrokAPI.py --list                         - List restore points
  python GrokAPI.py --search <keyword>             - Search history
  python GrokAPI.py --load [timestamp]             - Load history (all or specific)
  python GrokAPI.py --context <prompt>             - Call API with auto-context
  python GrokAPI.py <prompt> <model> <temp>        - Call with custom model/temp""",
            "response": ""
        }))
        sys.exit(1)
    
    command = sys.argv[1]
    
    # History management commands
    if command == "--list":
        result = list_restore_points()
        print(json.dumps(result, ensure_ascii=True, indent=2))
    
    elif command == "--search":
        if len(sys.argv) < 3:
            print(json.dumps({"success": False, "error": "Usage: --search <keyword>"}))
            sys.exit(1)
        keyword = sys.argv[2]
        result = search_history(keyword)
        # Truncate for display
        if result['success']:
            for m in result['matches']:
                m['content_snippet'] = m['content_snippet'][:200] + "..."
        print(json.dumps(result, ensure_ascii=True, indent=2))
    
    elif command == "--load":
        timestamp = sys.argv[2] if len(sys.argv) > 2 else None
        result = load_history(timestamp)
        # For display, just show structure not full content
        if result['success'] and result['count'] > 1:
            summary = {
                "success": True,
                "count": result['count'],
                "timestamps": result.get('timestamps', result.get('matched_timestamps', [])),
                "hint": "Use --load <timestamp> to view specific conversation"
            }
            print(json.dumps(summary, ensure_ascii=True, indent=2))
        else:
            print(json.dumps(result, ensure_ascii=True, indent=2))
    
    elif command == "--context":
        if len(sys.argv) < 3:
            print(json.dumps({"success": False, "error": "Usage: --context <prompt>"}))
            sys.exit(1)
        prompt = sys.argv[2]
        model = sys.argv[3] if len(sys.argv) > 3 else "grok-4-1-fast-reasoning"
        result = call_grok_api_with_context(prompt, model=model)
        print(json.dumps(result, ensure_ascii=True))
    
    else:
        # Standard API call
        prompt = sys.argv[1]
        model = sys.argv[2] if len(sys.argv) > 2 else "grok-4-1-fast-reasoning"
        temperature = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
        
        result = call_grok_api(prompt, model, temperature)
        # Use ensure_ascii=True to avoid UnicodeEncodeError on Windows cp1252 consoles
        # C++ QJsonDocument can handle escaped Unicode characters (\uXXXX)
        print(json.dumps(result, ensure_ascii=True))
