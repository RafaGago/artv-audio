#pragma once

#include <memory>

#include <JuceHeader.h>

extern juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter();

namespace artv {

class run_all_fx : public juce::ConsoleApplication::Command {
public:
  run_all_fx()
  {
    commandOption       = "--test-fx";
    argumentDescription = "";
    shortDescription    = "Runs all fx.";
    longDescription     = "";
    command = [this] (juce::ArgumentList const& args) { this->run (args); };
  }

private:
  void run (juce::ArgumentList const& args)
  {
    auto dsp = std::make_unique<juce::AudioProcessor> (createPluginFilter());
  }
};

} // namespace artv
