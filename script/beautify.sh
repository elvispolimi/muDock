#!/bin/bash

# check if git is available
git --version > /dev/null
rc=$?; if [[ $rc != 0 ]]; then echo "Please, install git" && exit $rc; fi

# check if clang-format is available
clang-format --version > /dev/null
rc=$?; if [[ $rc != 0 ]]; then echo "Please, install clang-format" && exit $rc; fi

# get the path to this script using the git repository facilities
if git rev-parse --git-dir > /dev/null 2>&1; then :
    PROJECT_ROOT=$(git rev-parse --show-toplevel)
else :
    echo "Error: you must execute this script in the mudock repository"
    exit -1
fi

# beautify the mudock library
clang-format -i -style=file $PROJECT_ROOT/mudock/include/mudock/*.hpp
clang-format -i -style=file $PROJECT_ROOT/mudock/include/mudock/*/*.hpp
clang-format -i -style=file $PROJECT_ROOT/mudock/src/*.cpp
clang-format -i -style=file $PROJECT_ROOT/mudock/src/*/*.cpp

# beautify the application
clang-format -i -style=file $PROJECT_ROOT/application/src/*.cpp
clang-format -i -style=file $PROJECT_ROOT/application/src/*.hpp
