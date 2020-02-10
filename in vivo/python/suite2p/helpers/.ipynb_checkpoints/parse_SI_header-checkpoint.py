import re
def parse_SI_header_level(header):
    rx_dict = {
        'stackNumSlices': re.compile(r'stackNumSlices = (?P<stackNumSlices>\d+)'),
    }
    for key, rx in rx_dict.items():
        match = rx.search(header)
        if match:
            level = int(match.group('stackNumSlices'))
            return level
    # if there are no matches
    return None

def parse_SI_header_zoom(header):
    rx_dict = {
        'scanZoomFactor': re.compile(r'scanZoomFactor = (?P<scanZoomFactor>\d+)'),
    }
    for key, rx in rx_dict.items():
        match = rx.search(header)
        if match:
            scanZoomFactor = int(match.group('scanZoomFactor'))
            return scanZoomFactor
    # if there are no matches
    return None

def parse_SI_header_FrameRate(header):
    rx_dict = {
        'scanFrameRate': re.compile(r'scanFrameRate = (?P<scanFrameRate>\d+)'),
    }
    for key, rx in rx_dict.items():
        match = rx.search(header)
        if match:
            scanFrameRate = int(match.group('scanFrameRate'))
            return scanFrameRate
    # if there are no matches
    return None

def parse_SI_header_Channels(header):
    rx_dict = {
        'channelsSave': re.compile(r'channelsSave = (?P<channelsSave>\d+)'),
    }
    for key, rx in rx_dict.items():
        match = rx.search(header)
        if match:
            scanFrameRate = int(match.group('channelsSave'))
            return scanFrameRate
    # if there are no matches
    return None

def parse_SI_header_Volumes(header):
    rx_dict = {
        'fastZNumVolumes': re.compile(r'fastZNumVolumes = (?P<fastZNumVolumes>\d+)'),
    }
    for key, rx in rx_dict.items():
        match = rx.search(header)
        if match:
            fastZNumVolumes = int(match.group('fastZNumVolumes'))
            return fastZNumVolumes
    # if there are no matches
    return None

def parse_SI_header_Frames(header):
    rx_dict = {
        'acqNumFrames': re.compile(r'acqNumFrames = (?P<acqNumFrames>\d+)'),
    }
    for key, rx in rx_dict.items():
        match = rx.search(header)
        if match:
            acqNumFrames = int(match.group('acqNumFrames'))
            return acqNumFrames
    # if there are no matches
    return None

def parse_SI_header_FramesPerFile(header):
    rx_dict = {
        'loggingFramesPerFile': re.compile(r'loggingFramesPerFile = (?P<loggingFramesPerFile>\d+)'),
    }
    for key, rx in rx_dict.items():
        match = rx.search(header)
        if match:
            loggingFramesPerFile = int(match.group('loggingFramesPerFile'))
            return loggingFramesPerFile
    # if there are no matches
    return None

