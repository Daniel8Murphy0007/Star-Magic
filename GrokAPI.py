#!/usr/bin/env python3
"""
GrokAPI.py - xAI Grok API Wrapper for SuperGrok4
Handles HTTPS requests to xAI API using Python requests library
Called by Source2.exe SuperGrok4Widget when Qt SSL is unavailable
"""

import sys
import os
import json
import requests

# Fix Windows console encoding issues (cp1252 -> UTF-8)
# Prevents UnicodeEncodeError when printing API responses with special characters
if sys.stdout.encoding != 'utf-8':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except AttributeError:
        # Python < 3.7 fallback
        import codecs
        sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')

def call_grok_api(prompt, model="grok-4-1-fast-reasoning", temperature=0.3):
    """
    Call xAI Grok API with given prompt
    
    Args:
        prompt: User question/prompt
        model: Grok model name (grok-4-1-fast-reasoning, grok-4-1-vision, grok-4-1)
        temperature: Response randomness (0.0-1.0)
    
    Returns:
        dict: {"success": bool, "response": str, "error": str}
    """
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

if __name__ == "__main__":
    # Command line usage: python GrokAPI.py <prompt> [model] [temperature]
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False,
            "error": "Usage: python GrokAPI.py <prompt> [model] [temperature]",
            "response": ""
        }))
        sys.exit(1)
    
    prompt = sys.argv[1]
    model = sys.argv[2] if len(sys.argv) > 2 else "grok-4-1-fast-reasoning"
    temperature = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
    
    result = call_grok_api(prompt, model, temperature)
    # Use ensure_ascii=True to avoid UnicodeEncodeError on Windows cp1252 consoles
    # C++ QJsonDocument can handle escaped Unicode characters (\uXXXX)
    print(json.dumps(result, ensure_ascii=True))
