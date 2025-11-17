# CMake generated Testfile for 
# Source directory: C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
# Build directory: C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test(CoreIntegrationTest "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/Debug/test_core_integration.exe")
  set_tests_properties(CoreIntegrationTest PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;292;add_test;C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test(CoreIntegrationTest "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/Release/test_core_integration.exe")
  set_tests_properties(CoreIntegrationTest PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;292;add_test;C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  add_test(CoreIntegrationTest "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/MinSizeRel/test_core_integration.exe")
  set_tests_properties(CoreIntegrationTest PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;292;add_test;C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test(CoreIntegrationTest "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/RelWithDebInfo/test_core_integration.exe")
  set_tests_properties(CoreIntegrationTest PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;292;add_test;C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/CMakeLists.txt;0;")
else()
  add_test(CoreIntegrationTest NOT_AVAILABLE)
endif()
