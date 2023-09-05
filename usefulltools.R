# In der Datei .Rbuildignore kann man Dateien (z.B. diese hier) und Ordner
# auflisten, die beim Zusammenpacken und Installieren des Paketes nicht
# beachtet werden sollen

## Folgendes erzeugt aus den Rogyxen-Kommentaren in den Funktionen automatisch
## die Hilfeseiten (im Ornder "man") und aktualisiert idR auch die NAMESPACE-Datei
## (die angibt, welche Funktionen innerhalb deines Paketes z.B. von anderen
## Paketen geladen werden und welche exportiert werden sollen)
devtools::document()

## Die Funktion läd alle Funktionen, die in deinem Paket drinne sind
## (so dass du schauen kannst, wie die Funktionen unter der Haube arbeiten können)
## (ruft glaube ich automatisch ein document() auf)
devtools::load_all()

#  Rhub hat viele Möglichkeiten das Paket zu testen. Das ist glaube ich aber
#  nur wichtig, wenn man Source-Code verwendet, also C++ oder so
rhub::check....

### Wichtige RStudio-Funktionen im "Build"-Panel
# Install: Guess what, installiert das Paket wie jedes andere
# Test: Ist was für später
# Check: Schaut ob das Paket in Ordnung ist, d.h. ob alle Funktionen und Variablen,
# die du innerhalb des Pakets verwendest auch definiert sind, ob du alle Argumente
# und deren Verwendung dokumentiert hast usw. gibt dir dann übersichtlich
# NOTES und ERRORs aus.
rhub::check_for_cran()

