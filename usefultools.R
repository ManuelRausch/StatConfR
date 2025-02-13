# In der Datei .Rbuildignore kann man Dateien (z.B. diese hier) und Ordner
# auflisten, die beim Zusammenpacken und Installieren des Paketes nicht
# beachtet werden sollen

## devtools::document() erzeugt aus den Rogyxen-Kommentaren in den Funktionen automatisch
## die Hilfeseiten (im Ornder "man") und aktualisiert idR auch die NAMESPACE-Datei
## (die angibt, welche Funktionen innerhalb deines Paketes z.B. von anderen
## Paketen geladen werden und welche exportiert werden sollen)
devtools::document()

#
devtools::build_manual()

urlchecker::url_check()

system("R CMD Rd2pdf C:/Users/PPA714/KU2/Projekte/StatConfR/StatConfR")

# this is the most important check! If you get errors here, it won't work with cran!
# it is much more critical than regular check!
system("R CMD check --as-cran ../statConfR_0.2.0.tar.gz")

# check on Rhub
rhub::rhub_setup()
rhub::rhub_doctor()
