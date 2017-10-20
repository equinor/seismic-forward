#include "nrlib/iotools/logkit.hpp"

#include "tasklist.hpp"

#include <iostream>

void TaskList::ViewAllTasks(void)
{
  NRLib::LogKit::WriteHeader("Suggested tasks");

  if (task_.size() > 0) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");
    for (size_t i=0 ; i < task_.size() ; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "%2d. %s\n", (i + 1), task_[i].c_str());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nNo tasks suggested.\n");
  }
}

std::vector<std::string> TaskList::task_(0);
