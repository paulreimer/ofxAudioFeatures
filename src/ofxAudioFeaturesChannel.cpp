/*
 * Copyright Paul Reimer, 2012
 *
 * This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
 * To view a copy of this license, visit
 * http://creativecommons.org/licenses/by-nc/3.0/
 * or send a letter to
 * Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
 */

#include "ofxAudioFeaturesChannel.h"

struct _aubio_onset_t
{
  aubio_pvoc_t * pv;        // phase vocoder
  aubio_specdesc_t * od;    // onset detection
  aubio_peakpicker_t * pp;  // peak picker
  cvec_t * fftgrain;        // phase vocoder output
  fvec_t * of;              // onset detection function
  smpl_t threshold;         // onset peak picking threshold
  smpl_t silence;           // silence threhsold
  uint_t minioi;            // minimum inter onset interval
  fvec_t * wasonset;        // number of frames since last onset
  uint_t samplerate;        // sampling rate of the input signal
  uint_t hop_size;          // number of samples between two runs
};

//--------------------------------------------------------------
ofxAudioFeaturesChannel::ofxAudioFeaturesChannel()
: gain(1.0)
, attack(0.3)
, decay(0.4)
, sampleRate(0.)
, hopSize(0)
, spectrumSize(0)
, bufferSize(0)
, currentHopBuffer(NULL)
, fftComplexOutputBuffer(NULL)
, fftProcessor(NULL)
, onsetOutputBuffer(NULL)
, onsetProcessor(NULL)
, pitchOutputBuffer(NULL)
, pitchProcessor(NULL)
, calibrateMic(false)
, calibratedMic(false)
{}

//--------------------------------------------------------------
ofxAudioFeaturesChannel::~ofxAudioFeaturesChannel()
{
  destroy();
}

//--------------------------------------------------------------
void
ofxAudioFeaturesChannel::setup(size_t _bufferSize, size_t _hopSize, float _sampleRate)
{
  bufferSize    = _bufferSize;
  hopSize       = _hopSize;
  sampleRate    = _sampleRate;
  spectrumSize  = _bufferSize / 2;

  fftComplexOutputBuffer  = new_cvec(bufferSize);
  fftProcessor            = new_aubio_fft(bufferSize);

  currentHopBuffer        = new_fvec(hopSize);
  onsetOutputBuffer       = new_fvec(1);
  onsetProcessor          = new_aubio_onset("specflux", bufferSize, hopSize, sampleRate);
  
  pitchOutputBuffer       = new_fvec(1);
  pitchProcessor          = new_aubio_pitch("yin", bufferSize, hopSize, sampleRate);
  aubio_pitch_set_unit(pitchProcessor, "bin");

	inputBuffer.resize(bufferSize);
/*
  spectralFeatureProcessor["energy"]    = new_aubio_specdesc("energy",    hopSize);
  spectralFeatureProcessor["hfc"]       = new_aubio_specdesc("hfc",       hopSize);
//  spectralFeatureProcessor["complex"]   = new_aubio_specdesc("complex",   hopSize);
//  spectralFeatureProcessor["phase"]     = new_aubio_specdesc("phase",     hopSize);
//  spectralFeatureProcessor["specdiff"]  = new_aubio_specdesc("specdiff",  hopSize);
//  spectralFeatureProcessor["kl"]        = new_aubio_specdesc("kl",        hopSize);
//  spectralFeatureProcessor["mkl"]       = new_aubio_specdesc("mkl",       hopSize);
//  spectralFeatureProcessor["specflux"]  = new_aubio_specdesc("specflux",  hopSize);
  spectralFeatureProcessor["centroid"]  = new_aubio_specdesc("centroid",  hopSize);
  spectralFeatureProcessor["spread"]    = new_aubio_specdesc("spread",    hopSize);
  spectralFeatureProcessor["skewness"]  = new_aubio_specdesc("skewness",  hopSize);
  spectralFeatureProcessor["kurtosis"]  = new_aubio_specdesc("kurtosis",  hopSize);
  spectralFeatureProcessor["slope"]     = new_aubio_specdesc("slope",     hopSize);
  spectralFeatureProcessor["decrease"]  = new_aubio_specdesc("decrease",  hopSize);
  spectralFeatureProcessor["rolloff"]   = new_aubio_specdesc("rolloff",   hopSize);
*/
//  for (std::map<std::string, aubio_specdesc_t*>::iterator spectralFeatureProcessorIter = spectralFeatureProcessor.begin();
//       spectralFeatureProcessorIter != spectralFeatureProcessor.end(); ++spectralFeatureProcessorIter)
  for (size_t i=0; i<usedFeatures.size(); ++i)
  {
    const std::string& featureName(usedFeatures[i]);
    std::vector<char> featureNameCopy(featureName.c_str(),
                                      featureName.c_str() + featureName.size() + 1);
    spectralFeatureProcessor.insert(make_pair(featureName,
                                              new_aubio_specdesc(&featureNameCopy[0], hopSize)));
    spectralFeatureOutputBuffer.insert(make_pair(featureName, new_fvec(1)));
  }

  spectrum.resize(spectrumSize);
	phase.resize(spectrumSize);
	binsScale.resize(spectrumSize);
  
	smoothedSpectrum.resize(spectrumSize, 1.);
	smoothedPhase.resize(spectrumSize, 1.);

	calibrateMic = false;
  
}

//--------------------------------------------------------------
void
ofxAudioFeaturesChannel::destroy()
{
  if (currentHopBuffer)
  {
    del_fvec(currentHopBuffer);
    currentHopBuffer = NULL;
  }

  if (fftProcessor)
  {
    del_aubio_fft(fftProcessor);
    fftProcessor = NULL;
  }

  if (fftComplexOutputBuffer)
  {
    del_cvec(fftComplexOutputBuffer);
    fftComplexOutputBuffer = NULL;
  }

  if (onsetProcessor)
  {
    del_aubio_onset(onsetProcessor);
    onsetProcessor = NULL;
  }

  if (onsetOutputBuffer)
  {
    del_fvec(onsetOutputBuffer);
    onsetOutputBuffer = NULL;
  }

  if (pitchProcessor)
  {
    del_aubio_pitch(pitchProcessor);
    pitchProcessor = NULL;
  }

  if (pitchOutputBuffer)
  {
    del_fvec(pitchOutputBuffer);
    pitchOutputBuffer = NULL;
  }

  for (std::map<std::string, aubio_specdesc_t*>::iterator spectralFeatureProcessorIter = spectralFeatureProcessor.begin();
       spectralFeatureProcessorIter != spectralFeatureProcessor.end(); ++spectralFeatureProcessorIter)
    del_aubio_specdesc(spectralFeatureProcessorIter->second);
  spectralFeatureProcessor.clear();
  
  for (std::map<std::string, fvec_t*>::iterator spectralFeatureOutputBufferIter = spectralFeatureOutputBuffer.begin();
       spectralFeatureOutputBufferIter != spectralFeatureOutputBuffer.end(); ++spectralFeatureOutputBufferIter)
    del_fvec(spectralFeatureOutputBufferIter->second);
  spectralFeatureOutputBuffer.clear();
}

#define lerp(start, stop, amt)  (start + (stop-start) * amt)
//--------------------------------------------------------------
void
ofxAudioFeaturesChannel::process(const float now)
{
  // input new hop
	for (unsigned int i=0; i<hopSize; ++i)
  {
		currentHopBuffer->data[i] = inputBuffer[i] * gain;
  }

  // process hop
  aubio_onset_do(onsetProcessor, currentHopBuffer, onsetOutputBuffer);
  aubio_pitch_do(pitchProcessor, currentHopBuffer, pitchOutputBuffer);

  for (std::map<std::string, aubio_specdesc_t*>::iterator spectralFeatureProcessorIter = spectralFeatureProcessor.begin();
       spectralFeatureProcessorIter != spectralFeatureProcessor.end(); ++spectralFeatureProcessorIter)
  {
    aubio_specdesc_do(spectralFeatureProcessorIter->second,
                      onsetProcessor->fftgrain, // should be pvoc output
                      spectralFeatureOutputBuffer[spectralFeatureProcessorIter->first]);

    spectralFeatures[spectralFeatureProcessorIter->first] = spectralFeatureOutputBuffer[spectralFeatureProcessorIter->first]->data[0];
  }

  // pitch extraction (per-hop)
  float currentPitch = pitchOutputBuffer->data[0];
  if (currentPitch != 0)
    pitch = lerp(pitch, currentPitch, 0.7);

  // onset extraction (per-hop)
  bool onsetDetected = (onsetOutputBuffer->data[0] > FLT_EPSILON);
  if (onsetDetected)
  {
    onsets.push_back(std::make_pair(now, pitch));
  }
  else if (!onsets.empty())
  {
    onsets.rbegin()->second = pitch;
  }

  // steal fft from onset detector's phase vocoder
  for (unsigned int i=0; i<spectrum.size(); ++i)
  {
    spectrum[i] = onsetProcessor->fftgrain->norm[i];
    phase[i] = onsetProcessor->fftgrain->phas[i];
  }

  // spectrum smoothing
  for (unsigned int i=0; i<spectrum.size(); ++i)
  {
    if (spectrum[i] > smoothedSpectrum[i])
      smoothedSpectrum[i] = lerp(spectrum[i], smoothedSpectrum[i], attack);
    else if (spectrum[i] < smoothedSpectrum[i])
      smoothedSpectrum[i] = lerp(spectrum[i], smoothedSpectrum[i], decay);
  }
  
  // spectrum magnitude calibration
  if (calibrateMic)
  {
    for (unsigned int i=0; i<smoothedSpectrum.size(); ++i)
    {
      if (binsScale[i] > FLT_MAX || binsScale[i] <= 0)
        binsScale[i] = 1.;
      else if (smoothedSpectrum[i] > binsScale[i])
        binsScale[i] = smoothedSpectrum[i];
    }

    calibratedMic = true;
  }
  else if (calibratedMic)
  {
    // scaling by calibrated values
    for (unsigned int i=0; i<smoothedSpectrum.size(); ++i)
      if (binsScale[i] > 0)
        smoothedSpectrum[i] /= binsScale[i];
  }
}
#undef lerp
