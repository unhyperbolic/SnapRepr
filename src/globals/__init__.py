
_defaultSettings = {}
_settings = {}

def getSetting(key):
    global _settings
    global _defaultSettings
    if _settings.has_key(key):
        return _settings[key]
    return _defaultSettings[key]

def registerSetting(key, defaultVal):
    global _defaultSettings
    _defaultSettings[key] = defaultVal

def setSetting(key, val):
    global _settings
    _settings[key] = val
