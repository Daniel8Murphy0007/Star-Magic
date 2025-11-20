std::cout << "\n=== REGISTRY ANALYSIS VIA WOLFRAM ===\n";
if (!InitializeWolframKernel()) {
    std::cout << "Kernel not available\n";
    return;
}

// Send all 810 registrations to Wolfram for analysis
std::string registrations = WolframEvalToString("Length[ReadList[\"COMPLETE_REGISTRY_LIST.txt\", String]]");
std::cout << "Total registrations in file: " << registrations << "\n";

// Analyze patterns
std::string analysis = WolframEvalToString(R"WL(
data = ReadList["COMPLETE_REGISTRY_LIST.txt", String];
{
  "Total" -> Length[data],
  "Unique" -> Length[DeleteDuplicates[data]],
  "SGR1745" -> Count[data, _?(StringContainsQ[#, "SGR1745"]&)],
  "SgrA" -> Count[data, _?(StringContainsQ[#, "SgrA"]&)],
  "Source" -> Count[data, _?(StringContainsQ[#, "Source"]&)],
  "Magnetar" -> Count[data, _?(StringContainsQ[#, "Magnetar"]&)],
  "NGC" -> Count[data, _?(StringContainsQ[#, "NGC"]&)]
}
)WL");

std::cout << "Wolfram Analysis:\n" << analysis << "\n";
