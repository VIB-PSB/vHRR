import configparser
from pathlib import Path

def read_config():
    config = configparser.ConfigParser()
    try:
        config.read('{}/.gbaconfig'.format(Path.home()))
    except:
        raise FileNotFoundError('[ERROR] GBA config file missing: run scripts/gba_config.py (usage: python3 gba_config.py [options])')
    return config

def get_obo():
    config = read_config()
    return config['GO']['go_obo']