# API Keys Setup Instructions

## Quick Setup

1. **Copy the example environment file:**
   ```powershell
   Copy-Item .env.example .env
   ```

2. **Edit the `.env` file and add your xAI API key:**
   ```
   XAI_API_KEY=your_xai_key_here
   ```

3. **Load environment variables (PowerShell):**
   ```powershell
   # Manual load for current session (replace with your actual key)
   $env:XAI_API_KEY = "your_xai_key_here"
   
   # Or load all from .env file
   Get-Content .env | ForEach-Object {
       if ($_ -match '^([^=]+)=(.*)$') {
           [Environment]::SetEnvironmentVariable($matches[1], $matches[2], 'Process')
       }
   }
   ```

4. **Verify the key is loaded:**
   ```powershell
   echo $env:XAI_API_KEY
   ```

## Configured API Keys

### NASA APIs (Already in code as fallback)
- **NASA_API_KEY_1**: `PNJaNeFWqMb2g0CEQGqJePkndqYfKvBzq6XJqAwg`
- **NASA_API_KEY_2**: `FJnBo64nLFqExHwDchrcaf101D8wmGSm0cF27clz`
- **MAST_API_KEY**: `emXvt90Htf0U4RogKTB5lqSxClUeg2pvMQxvZciM`

### xAI (Grok) - For AI Fallback
- **XAI_API_KEY**: Set via environment variable (see `.env.example`)
  - Used when SIMBAD/NED don't have complete data
  - Provides intelligent estimates for missing parameters

### Optional APIs
- **OPENAI_API_KEY**: Alternative AI fallback
- **ANTHROPIC_API_KEY**: Alternative AI fallback (Claude)
- **WOLFRAM_APP_ID**: For symbolic computation and verification

## Python Usage

```python
from APIFetch import fetch

# The environment variable will be automatically loaded
result = fetch("Sagittarius A*")

# Grok will be used as fallback if SIMBAD/NED are incomplete
```

## Security Notes

- ✅ `.env` is in `.gitignore` and will NOT be committed
- ✅ API keys are loaded from environment variables
- ✅ `.env.example` contains the working xAI key for easy setup
- ⚠️ Never commit `.env` file with actual keys to public repositories
