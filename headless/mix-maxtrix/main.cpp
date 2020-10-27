#if 0 // went for smoke tests...

#include <JuceHeader.h>

#include "mix-maxtrix/run_all_fx.hpp"

using namespace juce;

String getVersion()
{
  return String (ProjectInfo::projectName) + " - " + ProjectInfo::versionString;
}

String getHelp()
{
  return "mix-maxtrix-4 smoke tester:";
}

int main (int argc, char const* argv[])
{
#if JUCE_MAC
  Process::setDockIconVisible (false); // hide dock icon
#endif
  ScopedJuceInitialiser_GUI scopedJuce; // creates MessageManager

  ConsoleApplication app;
  app.addVersionCommand ("--version", getVersion());
  app.addHelpCommand ("--help|-h", getHelp(), true);

  return app.findAndRunCommand (argc, argv);
}

#endif
