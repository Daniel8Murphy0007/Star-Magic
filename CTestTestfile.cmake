# CMake generated Testfile for 
# Source directory: C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
# Build directory: C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test(SOURCE4_DualMethodValidation "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/Debug/test_source4_validation.exe")
  set_tests_properties(SOURCE4_DualMethodValidation PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;1047;add_test;C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test(SOURCE4_DualMethodValidation "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/Release/test_source4_validation.exe")
  set_tests_properties(SOURCE4_DualMethodValidation PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;1047;add_test;C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  add_test(SOURCE4_DualMethodValidation "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/MinSizeRel/test_source4_validation.exe")
  set_tests_properties(SOURCE4_DualMethodValidation PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;1047;add_test;C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test(SOURCE4_DualMethodValidation "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/RelWithDebInfo/test_source4_validation.exe")
  set_tests_properties(SOURCE4_DualMethodValidation PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;1047;add_test;C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;0;")
else()
  add_test(SOURCE4_DualMethodValidation NOT_AVAILABLE)
endif()
