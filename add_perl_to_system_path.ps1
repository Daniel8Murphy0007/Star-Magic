# Run as Administrator to add Perl to System PATH
$currentPath = [Environment]::GetEnvironmentVariable("PATH", "Machine")
$perlPath = "C:\Strawberry\perl\bin"

if ($currentPath -notlike "*$perlPath*") {
    [Environment]::SetEnvironmentVariable("PATH", "$perlPath;$currentPath", "Machine")
    Write-Output "SUCCESS: Added $perlPath to System PATH"
    Write-Output "Please close and reopen terminals for changes to take effect"
}
else {
    Write-Output "Perl already in System PATH"
}
