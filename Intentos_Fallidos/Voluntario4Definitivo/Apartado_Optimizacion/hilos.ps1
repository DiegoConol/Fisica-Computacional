# Compila y ejecuta optimizadoomp.c con diferentes números de hilos usando OpenMP

$source = "optimizadoomp.c"
$out = "optimizadoomp.exe"

# Compilar con OpenMP
Write-Host "Compilando con OpenMP..."
gcc $source -o $out -fopenmp

# Lista de hilos a probar
$hilos = @(1, 2, 4, 8, 12)

foreach ($n in $hilos) {
    Write-Host "`nEjecutando con $n hilo(s)..."
    $env:OMP_NUM_THREADS = $n
    Measure-Command { .\optimizadoomp.exe } | Select-Object TotalSeconds
}