#pragma once

#import <AudioToolbox/AudioServices.h>

class ofxAudioDeviceControl
{
public:
//  static size_t getOutputChannels();
  static float getVolume();
  static void setVolume(float newVolume);
  static AudioDeviceID defaultOutputDeviceID();
  
  static size_t getInputChannels();
  static float getInputVolume();
  static void setInputVolume(float newVolume);
  static AudioDeviceID defaultInputDeviceID();
};

