# Compila con diferentes optimizaciones y OpenMP, luego ejecuta cada uno y mide el tiempo

$source = "sinoptimizar.c"

# Compilaciones
Write-Host "Compilando con OpenMP..."
gcc $source -o programaOpen -fopenmp

Write-Host "Compilando con OpenMP y -O1..."
gcc $source -o programaO1 -fopenmp -O1

Write-Host "Compilando con OpenMP y -O2..."
gcc $source -o programaO2 -fopenmp -O2

Write-Host "Compilando con OpenMP y -O3..."
gcc $source -o programaO3 -fopenmp -O3

Write-Host "Compilando con OpenMP y -Ofast..."
gcc $source -o programaOfast -fopenmp -Ofast

# Ejecuciones y medición de tiempo
Write-Host "`nEjecutando programaOpen..."
Measure-Command { .\programaOpen.exe } | Select-Object TotalSeconds

Write-Host "`nEjecutando programaO1..."
Measure-Command { .\programaO1.exe } | Select-Object TotalSeconds

Write-Host "`nEjecutando programaO2..."
Measure-Command { .\programaO2.exe } | Select-Object TotalSeconds

Write-Host "`nEjecutando programaO3..."
Measure-Command { .\programaO3.exe } | Select-Object TotalSeconds

Write-Host "`nEjecutando programaOfast..."
Measure-Command { .\programaOfast.exe } | Select-Object TotalSeconds