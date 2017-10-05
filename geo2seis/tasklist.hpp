#ifndef TASKLIST_H
#define TASKLIST_H

#include <vector>
#include <string>

class TaskList
{
public:
  static void AddTask(const std::string & task) { task_.push_back(task) ;}
  static void ViewAllTasks(void);

private:
  static std::vector<std::string> task_;
};

#define TASKLIST_END
#elif !defined TASKLIST_END
#error tasklist.h is part of a cyclic dependency structure
#endif
