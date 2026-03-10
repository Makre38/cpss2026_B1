set terminal gif animate delay 8 size 640,640
set output "output/ising.gif"

unset key
unset colorbox
unset xtics
unset ytics
set size square
set autoscale fix
set lmargin at screen 0.08
set rmargin at screen 0.92
set bmargin at screen 0.08
set tmargin at screen 0.90
set cbrange [-1:1]
set palette defined (-1 "#4169E1", 0 "#FFFFFF", 1 "#FF6347")

files = system("find output/snapshots -name 'sweep_*.dat' | sort")
nfiles = words(files)

do for [i=1:nfiles] {
    file = word(files, i)
    set title sprintf("2D Ising model: sweep %d", i - 1)
    plot file matrix with image pixels
}
