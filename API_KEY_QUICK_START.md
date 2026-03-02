# SuperGrok4 API Key Quick Start Guide

## ⚡ 30-Second Setup

1. **Get API Key** → Visit [https://x.ai/api](https://x.ai/api) → Sign up → Generate key
2. **Open source2.exe** → Go to Tab 7 (SuperGrok4)
3. **Click "Configure API Key"** → Paste key in text field
4. **Click "💾 Save API Key to Config"** → Done!
5. **Next time**: Key loads automatically 🎉

That's it! No environment variables, no app restart needed.

---

## 📍 Where Button Is Located

**source2.cpp Tab 7 - SuperGrok4 Physics Assistant**
- Look for blue button labeled: **"🔧 Configure API Key"**
- It's in the top toolbar with other settings

---

## 🔑 What Your API Key Looks Like

✅ **Valid Format**: `xai-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx` (starts with `xai-`)

❌ **Invalid Format**: Anything else

---

## 💾 Where Is My Key Saved?

**File**: `grok_api_config.json` in your Star-Magic project folder

**Format**: Plain JSON (human-readable)

**Located**: `C:\Users\...\source\repos\Daniel8Murphy0007\Star-Magic\grok_api_config.json`

---

## ❓ Common Issues

| Problem | Solution |
|---------|----------|
| "Key didn't save" | Make sure you clicked "💾 Save API Key to Config" button |
| "Still shows ❌ after save" | Restart source2.exe or refresh dialog |
| "Key not working" | Verify it starts with `xai-` and is valid from x.ai/api |
| "Can't find the dialog" | Go to Tab 7, look for 🔧 icon in toolbar |
| "Dialog says 'APIKeyManager not found'" | APIKeyManager.py must be in same folder as source2.cpp |

---

## 🔄 Switching Keys

1. Click "Configure API Key"
2. Enter new key in text field
3. Click "💾 Save API Key to Config"
4. Old key replaced automatically

---

## 🗑️ Removing Saved Key

1. Click "Configure API Key"
2. Leave text field empty
3. Click "🗑️ Clear Saved Key"
4. Confirm deletion

---

## 🌐 Models You Can Use

All these models work with your API key:

- **grok-4-1-fast-reasoning** (default, fastest)
- **grok-4-1** (best quality)
- **grok-4-1-vision** (can analyze images)

---

## 💡 Pro Tips

1. **Keep Key Secret**: Don't share `grok_api_config.json` or commit to public repos
2. **Add to .gitignore**: `echo "grok_api_config.json" >> .gitignore`
3. **Environment Variable Fallback**: Still works if you prefer `$env:XAI_API_KEY = "..."`
4. **Multiple Computers**: Copy `grok_api_config.json` to the same folder on other machines
5. **Security**: For production, consider encrypting the config file

---

## 🆘 Need Help?

- **Grok API Website**: https://x.ai/api
- **API Key Issues**: Check your key at https://console.x.ai/api-keys
- **Application Issues**: See `API_KEY_MANAGEMENT_FEATURE.md` for full documentation

---

**Version**: 1.0 | **Last Updated**: March 2, 2026

