#!/usr/bin/env python3
"""
OllamaAPI.py - Ollama API Wrapper for CoAnQi_bot
Handles HTTP requests to local Ollama server (localhost:11434)
Called by Source2.exe OllamaCodeBotWidget
"""

import sys
import os
import json
import requests

# Fix Windows console encoding issues (cp1252 -> UTF-8)
if sys.stdout.encoding != 'utf-8':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except AttributeError:
        # Python < 3.7 fallback
        import codecs
        sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')

def call_ollama_api(prompt, model="llama3.2", system_prompt=None):
    """
    Call Ollama API with given prompt
    
    Args:
        prompt: User question/prompt
        model: Ollama model name (llama3.2, codellama, qwen2.5-coder, etc.)
        system_prompt: Optional system prompt to set bot behavior
    
    Returns:
        dict: {"success": bool, "response": str, "error": str}
    """
    
    url = "http://localhost:11434/api/generate"
    
    if system_prompt is None:
        system_prompt = (
            "You are CoAnQi_bot, an expert physics and code assistant for the UQFF "
            "(Unified Quantum Field Framework) project. You have deep knowledge of "
            "astrophysics, quantum mechanics, UQFF equations, and scientific computing. "
            "Help with code generation, physics equations, and scientific computing. "
            "Provide detailed explanations with equations and code examples when relevant. "
            "Be precise and comprehensive."
        )
    
    payload = {
        "model": model,
        "prompt": prompt,
        "system": system_prompt,
        "stream": False
    }
    
    try:
        response = requests.post(url, json=payload, timeout=120)
        
        if response.status_code == 200:
            data = response.json()
            bot_response = data.get("response", "")
            
            if not bot_response:
                return {
                    "success": False,
                    "error": "Empty response from Ollama",
                    "response": ""
                }
            
            return {
                "success": True,
                "response": bot_response,
                "error": ""
            }
        
        elif response.status_code == 404:
            # Model not found
            error_data = response.json()
            error_msg = error_data.get("error", f"Model '{model}' not found")
            return {
                "success": False,
                "error": f"{error_msg}. Install with: ollama pull {model}",
                "response": ""
            }
        
        else:
            # Other HTTP errors
            try:
                error_data = response.json()
                error_msg = error_data.get("error", f"HTTP {response.status_code}")
            except:
                error_msg = f"HTTP {response.status_code}: {response.text[:200]}"
            
            return {
                "success": False,
                "error": error_msg,
                "response": ""
            }
    
    except requests.exceptions.ConnectionError:
        return {
            "success": False,
            "error": "Cannot connect to Ollama. Ensure 'ollama serve' is running on localhost:11434",
            "response": ""
        }
    
    except requests.exceptions.Timeout:
        return {
            "success": False,
            "error": "Request timeout (120 seconds) - model may be loading or prompt too complex",
            "response": ""
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": f"Unexpected error: {str(e)[:200]}",
            "response": ""
        }

if __name__ == "__main__":
    # Command line usage: python OllamaAPI.py <prompt> [model] [system_prompt]
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False,
            "error": "Usage: python OllamaAPI.py <prompt> [model] [system_prompt]",
            "response": ""
        }))
        sys.exit(1)
    
    prompt = sys.argv[1]
    model = sys.argv[2] if len(sys.argv) > 2 else "llama3.2"
    system_prompt = sys.argv[3] if len(sys.argv) > 3 else None
    
    result = call_ollama_api(prompt, model, system_prompt)
    # Use ensure_ascii=True to avoid UnicodeEncodeError on Windows cp1252 consoles
    # C++ QJsonDocument can handle escaped Unicode characters (\uXXXX)
    print(json.dumps(result, ensure_ascii=True))
