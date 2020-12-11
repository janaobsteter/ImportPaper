## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Seštejemo števili 3 + 5
3 + 5


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Spremenljivki stevilo pripišemo vrednosti 3
stevilo <- 3


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Nastavimo stevilo
# stevilo <- 3


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Kot primer naložimo paket 'base', ki je kot privzet vključen v R
library(base)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo spremenljivko "stevilo" z vrednostjo 10
stevilo <- 10
# Izpišimo "stevilo"
stevilo
# Preverimo tip spremenljivke
class(stevilo)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo spremenljivko "stevilo" z vrednostjo 10 kot "integer"
stevilo <- as.integer(10)
# Preverimo tip spremenljivke
class(stevilo)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Direktno ustvarimo "integer" spremenljivko "stevilo" z vrednostjo 10
stevilo <- 10L
class(stevilo)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo spremenljivko "realnoStevilo" z vrednostjo 3.52
realnoStevilo <- 3.52
# Izpišimo "realnoStevilo"
realnoStevilo
# Preverimo tip spremenljivke
class(realnoStevilo)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo spremenljivko "beseda" z vrednostjo "jagoda"
beseda <- "jagoda"
# Izpišimo spremenljivko "beseda"
beseda
# Preverimo tip spremenljivke
class(beseda)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo spremenljivko "logicnaVrednost" z vrednostjo TRUE
logicnaVrednost <- TRUE
# Izpišimo "logicnaVrednost"
logicnaVrednost
# Preverimo tip spremenljivke
class(logicnaVrednost)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Sesštejmo vrednosti TRUE in TRUE
TRUE + TRUE


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo spremenljivko "manjkaVrednost" z manjkajočo vrednostjo
manjkaVrednost <- NA
# Izpišimo "manjkaVrednost"
manjkaVrednost
# Preverimo tip spremenljivke
class(manjkaVrednost)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo prazen element "prazenElement"
prazenElement <- NULL
# Izpišimo "prazenElement"
prazenElement
# Preverimo tip spremenljivke
class(prazenElement)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo številski vektor "stevila" z vrednostmi 1, 10 in 100
stevila <- c(1, 10, 100)
# Izpišimo vektor
stevila


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo številski vektor "meseci" s številkami mesecev 1 do 12
meseci <- 1:12
# Izpišimo vektor
meseci


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo vektor znakov "letniCasi" z letnimi časi
letniCasi <- c("pomlad", "poletje", "jesen", "zima")
# Izpišimo vektor
letniCasi



## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo drugi element vektorja "letniCasi"
letniCasi[2]


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo faktor "letniCasi_f" s ponovljenimi vrednostmi za letne čase 
# c("pomlad", "pomlad", "poletje", "poletje", "jesen", "jesen", "zime", "zima")
letniCasi_f <- factor(c("pomlad", "pomlad", "poletje", "poletje", "jesen", "jesen", "zime", "zima"))
# Izpišimo "letniCasi_f"
letniCasi_f


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Poskušajmo dodati novo raven v faktor
letniCasi_f[9] <- "zjutraj"
# Ta operacija je nedovoljena, saj "zjutraj" ni ena izmed ravni faktorja (ne moremo ustvariti novih)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo seznam "leto", ki vključuje vektorje 
# "meseci" in "letniCasi" iz prejsnjega koraka in število 2020
leto <- list(meseci, letniCasi, 2020)
# Izpišimo seznam "leto"
leto

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo dolžino seznama "leto""
length(leto)



## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberemo prvi element seznama "leto"
leto[[1]]

## -----------------------------------------------------------------------------------------------------------------------------------------------------

# Ker je prvi element vektor, izberimo še peti element znotraj tega
leto[[1]][5]


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo tabelo "povpTemp" z letnimi časi in povprečno temperature za vsakega izmed njih
# Predpostavimo, da je povprečna temperatura pomladi je 18, poleti 24, jeseni 13 in pozimi 5
povpTemp <- data.frame(LetniCas = letniCasi, temp = c(18, 24, 13, 5))
# Izpišimo tabelo
print(povpTemp)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izpišimo prvih par vrstic tabele "povpTemp"
head(povpTemp)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izpišimo prvi dve vrstici tabele "povpTemp"
head(povpTemp, n=2)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo število vrstic in stolpcev tabele "povpTemp"
nrow(povpTemp)
ncol(povpTemp)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo imena vrstic in stolpcev tabele "povpTemp"
rownames(povpTemp)
colnames(povpTemp)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo prvo vrednost v prvem stolpcu tabele "povpTemp"
povpTemp[1,1]

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo prvo vrstico tabele "povpTemp"
povpTemp[1,]

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo prvi stolpec tabele "povpTemp"
povpTemp[,1]


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo stolpec <LetniCasi> v tabeli "povpTemp"
povpTemp$LetniCas


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo matriko "vrednosti" z vrednostmi 1:20, štirimi vrsticami in petimi stolpci
vrednosti <- matrix(1:20, nrow=4, ncol=5)
# Izpišimo matriko
vrednosti


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo element v drugi vrstici tretjega stolpca matrike vrednosti
vrednosti[2, 3]


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo polje povpTempLeto s povprečnimi letnimi temperaturami v 
# letnih časih v dveh različnih letih.
# Povprečna temparatura je vektor vrednosti 18, 24, 13, 5, 15, 27, 12, 7; leti sta 2018 in 2019
# Dimenzije: x = letni čas (4), y = povprečna temperatura (1), z = leto (2)
povpTempLeto <- array(c(18, 24, 13, 5, 15, 27, 12, 7), 
                      dim=c(4, 1, 2), 
                      dimnames=list(letniCasi,
                                    "povpTemp",
                                    c(2018, 2019)))
# Izpišimo polje "povpTempLeto"
povpTempLeto


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo dimenzije polja "povpTempLeto"
dim(povpTempLeto)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo strukturo polja "povpTempLeto"
str(povpTempLeto)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo spremenljivki "a" z vrednostjo 2 in "b" z vrednostjo 11.2
a <- 2
b <- 11.2
# Ustvarimo spremenvljivko "vsota" kot vsoto a in b
(vsota <- a + b)
# Ustvarimo spremenvljivko "razlike" kot vsoto a in b
(razlika <- a - b)
# Ustvarimo spremenvljivko "zmnožek" kot vsoto a in b
(zmnožek <- a * b)
# Ustvarimo spremenvljivko "kvocient" kot vsoto a in b
(kvocient <- a / b)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo vektorja "a" z vrednostmi c(1, 2, 3) in "b" z vrednostmi c(3.12, 5.44, 10)
a <- c(1, 2, 3) 
b <- c(3.12, 5.44, 10)
# Ustvarimo spremenvljivko "vsota" kot vsoto a in b
(vsota <- a + b)
# Ustvarimo spremenvljivko "razlike" kot vsoto a in b
(razlika <- a - b)
# Ustvarimo spremenvljivko "zmnožek" kot vsoto a in b
(zmnožek <- a * b)
# Ustvarimo spremenvljivko "kvocient" kot vsoto a in b
(kvocient <- a / b)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo ostanek pri deljenju 10 / 3
10 %% 3
# Preverimo celoštevilski kvocient deljenja 10 / 3
10 %/% 3
# Ustvarimo matriko m1 z vrednostmi 1:9 in tremi vrstiami 
# ter matriko m2 z vrednostmi 9:18 in tremi vrsticami
m1 <- matrix(1:9, nrow=3)
m2 <- matrix(10:18, nrow=3)
m1 %*% m2


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo spremenljivko a z vrednotjo 5
a <- 5
# Preverimo, ali je 5 večje od 3
a > 3

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo, ali je 5 enako 3
a == 3

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo, ali 5 ni enako 3
a != 3

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo, ali vektor c(1, 2, 3) vsebuje 4
a %in% c(1,2,3)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Previmo, koliko vrednosti a (5) je v vektorju vrednosti = c(1, 3, 5, 4, 5, 7)
vrednosti = c(1, 3, 5, 4, 5, 7)
sum(vrednosti == a)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Združimo besedi Janez in Novak s presledkom
paste("Janez", "Novak", sep=" ")


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Združimo besedi Janez in Novak brez separatorja
paste0("Janez", "Novak")


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo vektor "stevila" z vrednostmi 1:10
stevila <- 1:10
# Seštejmo elemente vektorja s funkcijo sum()
(sum(stevila))

## -----------------------------------------------------------------------------------------------------------------------------------------------------
#Izračunajmo povprečje in varianco vektorja "stevila" s funkcijo mean() in var()
(mean(stevila))
(var(stevila))


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo dolžino vektorja stevila s funkcijo length()
(length(stevila))
# Preverimo minimalno in maksimalno vrednost vektorja "stevila" s funkcijama min() in max()
(min(stevila))
(max(stevila))


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo vektor vrednosti <- c(1, 2, 3, 4, NA)
vrednosti <- c(1, 2, 3, 4, NA)
# Seštejmo vektor vrednosti
sum(vrednosti)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo pomoč za funkcijo sum()
?sum()
# Ponovno sštejmo vektor vrednosti
# Ker želimo odstraniti manjkajočih vrednosti, nastavimo parameter na.rm na FALSE
(sum(vrednosti, na.rm = TRUE))


## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------
# Nastavimo delovni imenik na lokalno mapo z datotekami
workingDir <- "/home/jana/Documents/Delavnice/R/Material/"
knitr::opts_knit$set(root.dir = '/home/jana/Documents/Delavnice/R/Material/')
setwd(workingDir)



## -----------------------------------------------------------------------------------------------------------------------------------------------------
getwd()


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preberimo datoteko ZivljenjskaDoba_mali.csv
zivDoba <- read.csv("ZivljenjskaDoba.csv")


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo prvih par vrstic prebrane tabele s funkcijo head()
head(zivDoba)

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo zadnjih par vrstic prebrane tabele s funkcijo tail()
tail(zivDoba)

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo strukturo prebrane tabele s funkcijo str()
str(zivDoba)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Vidimo, da je leto kodirano kot številska spremenljivka (integer). 
# S funkcijo as.factor() jo spremenimo v faktor (kategorično spremenljivko).
zivDoba$Leto <- as.factor(zivDoba$Leto)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo število vrstic in število stolpcev
nrow(zivDoba)
ncol(zivDoba)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preimenujmo stolpec "Country" v "Drzava"
colnames(zivDoba)[1] <- "Drzava"


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo tretjo vrstico tabele zivDoba z indeksom
zivDoba[3,]

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo tretji stolpec tabele zivDoba z indeksom [vrstica, stolpec]
zivDoba[, 3]

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo tretji stolpcev tabele zivDoba z imenom (Status)
zivDoba$Status



## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo, koliko podatkov (vrstic) imamo po državi s funkcijo table()
table(zivDoba$Drzava)



## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo vsebnost stolpca ZivDoba s funkcijo summary
summary(zivDoba$ZivDoba)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo povprečje življenjske dobe preko vseh držav
mean(zivDoba$ZivDoba)
# Preverimo varianco življenjske dobe
var(zivDoba$ZivDoba)
# Preverimo minimalno umrljivost v podatkih
min(zivDoba$Umrljivost)
# Preverimo maksimalno umrljivost v podatkih
max(zivDoba$Umrljivost)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo vse vrstice tabele, v katerih je vrednost spremenljivke ZivDoba večja od 85
zivDoba[zivDoba$ZivDoba > 85,]

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Izberimo vse vrstice, v katerih je vrednost stolpca "Drzava" enaka "France"
zivDoba[zivDoba$Drzava == "France",]


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Odstranimo stolpec "Polio" z indeksom iz tabele zivDoba
# novo ustvarjeno tabelo poimenujemo zivDoba_nov
zivDoba_nov <- zivDoba[, -10]
# Preverimo imena stolpcev nove tabele zivDoba_nov s funkcijo colnames()
colnames(zivDoba_nov)



## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo vektor 'ospiceFR', v katerega zberemo podatke o letu in pojavnosti ošpic v Franciji
ospiceFR <- zivDoba[zivDoba$Drzava == "France", c("Leto", "Ospice")]
head(ospiceFR)

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preimenuj drugi stolpceviz 'Ospice' v 'ospiceFR'
colnames(ospiceFR)[2] <- 'ospiceFR'
head(ospiceFR)

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo vektor 'ospiceBR', v katerega zberemo podatke o letu in pojavnosti ošpic v Braziliji
ospiceBR <- zivDoba[zivDoba$Drzava == "Brazil", c("Leto", "Ospice")]
# Preimenuj drugi stolpceviz 'Ospice' v 'ospiceBR'
colnames(ospiceBR)[2] <- 'ospiceBR'


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo vektor 'ospiceAU' in' ospiceEG'
# Vanju zberemo podatke o pojavnosti ošpic v Avstriji in Egiptu
ospiceAU <- zivDoba$Ospice[zivDoba$Drzava == "Austria"]
ospiceEG <- zivDoba$Ospice[zivDoba$Drzava == "Egypt"]


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarimo nov stolpec "Celina" v tabeli zivDoba in mu pripišimo vrednost "Evropa"
zivDoba$Celina <- "Evropa"
# Preverimo tabelo
head(zivDoba)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Vse države v tabeli zivDob ne ležijo v Evropi, in sicer "Egypt", "Brazil", "Chile" in "Peru"
# Naprej v vsrticah z vrednostjo "Egypt" spremenimo vrednost v stolpcu "Celina" na "Afrika"
# 
zivDoba$Celina[zivDoba$Drzava == "Egypt"] <- "Afrika"

## -----------------------------------------------------------------------------------------------------------------------------------------------------
#"Brazil", "Chile" in "Peru" ležijo v Južni Ameriki
# Tem vrsticam spremenimo vrednost v stolpcu "Celina" na JuznaAmerika
zivDoba$Celina[zivDoba$Drzava %in% c("Brazil", "Peru", "Chile")] <- "JuznaAmerika"

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Prevrimo zastopanost vrednosti v stolpcu "Celina"
table(zivDoba$Celina)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo write.csv() zapišemo popravljeno tabelo zivDoba v datoteko 
# ZivljenjskaDoba_popravljena.csv
# Shranimo podatke brez navednic in brez imena vrstic
write.csv(zivDoba, "ZivljenjskaDoba_popravljena.csv", quote=FALSE, row.names=FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo cbind() združimo zgoraj ustvarjene vektorje 'ospiceFR', 'ospiceAU' in 'ospiceEG' 
# v novo tabelo 'ospice'
ospice <- cbind(ospiceAU, ospiceEG)
# Preverimo novonastalo tablelo
head(ospice)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo merge() združimo zgoraj ustvarjene vektorje 'ospiceFR' in 'ospiceBR' 
# na podlagi spremenljivke 'Leto' v novo tabelo 'ospiceMerged'
ospiceMerged <- merge(ospiceFR, ospiceBR, by="Leto")
# Preverimo novonastalo tablelo
head(ospiceMerged)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Recimo, da želimo imeti en stolpec z vsemi vrednostmi za hepatitisB, polio in HIV za vse države, 
# v drugem stolpcu pa bi imeli spremenljivko, za katero bolezen gre
# Takšno organizacijo podatkov ponavadi potrebujemo za grafično prikazovanje podatkov
# Najprej izberemo željene stolpce celotne tabele zivDoba - Drzava, Leto, HepatitisB, HIV in Polio
bolezni <- zivDoba[, c("Drzava", "Leto", "HepatitisB", "HIV", "Polio")]
# Preverimo novonastalo tabelo
head(bolezni)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo melt() iz paketa reshape preuredumo tabelo, 
# tako da sta naša ključa Leto in Drzava (stolpca bosta ostala nespremenjena)
# Novo nastalo tabelo poimenujemo bolezni_melt
library(reshape)
bolezni_melt <- melt(bolezni, id.vars = c("Leto", "Drzava"))
# Preveimo novonastalo tabelo s funkcijo head()
head(bolezni_melt)

## -----------------------------------------------------------------------------------------------------------------------------------------------------
#Preverimo vrednosti v stolpcu 'variable' tabele bolezni_melt
table(bolezni_melt$variable)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ponovimo vajo
# Sfunkcijo pivot_longer() iz paketa tidyr preuredimo tabelo bolezni (ustvarjena zgoraj), 
# tako da zberemo vse vredosti za stolpce HepatitisB, HIV in Polio v enem stolpcu
# Novo nastalo tabelo poimenujemo bolezni_longer
library(tidyr)
bolezni_longer <- pivot_longer(data = bolezni, cols = c(HepatitisB, HIV, Polio), names_to = "Bolezen")
# Preverimo novonastalo tabelo
head(bolezni_longer)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo pivot_wider() razširimo zgornjo tabelo bolezni_longer, 
# da bodo podatki za vsako državo v svojem stolpcu
bolezni_wider <- pivot_wider(data = bolezni_longer, id_cols = c(Leto, Bolezen), 
                             names_from=Drzava, values_from = value)
# Preverimo novonastalo tabelo
head(bolezni_wider)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Združimo stolpca Drzava in Leto tabele bolezni z ločilom "_" v nov stolpev DrzavaLeto
head(unite(bolezni, DrzavaLeto, Drzava, Leto, sep="_"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo expand() preverimo vse kombinacije spremenljivk Drzava in Leto 
# v tabelo bolezni in rezultat shranimo v bolezni_kombo
bolezni_kombo <- expand(bolezni, Drzava, Leto)
# Preverimo število kombinacij
nrow(bolezni_kombo)
# Prikažemo zadnjih par vrstic tabele
tail(bolezni_kombo)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Uporabimo funkcijo chop(), da vse vrednosti za bolezni po letih shranimo 
# v eno spremenljivko glede na Drzavo. Rezultat shranimo v tabelo bolezni_chop.
bolezni_chop <- chop(bolezni, c(Leto, HepatitisB, HIV, Polio))
# Preverimo rezultat
head(bolezni_chop, 3)

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo dimenzije vstopne in novonastale tabele
dim(bolezni)
dim(bolezni_chop)
# Novonastala tabelo je dimenzije


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Novonastala tabela vsebuje samo eno vrstico za vsako državo!
# Elementi tabele so seznami vektorjev
# Preverimo stolpec leto tabele bolezni_chop
head(bolezni_chop$Leto, 3)

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Preverimo stolpev Polio tabele bolezni_chop
head(bolezni_chop$Polio, 3)
# Elementi so seznami z vektorji za vsak nivo preostalih spremenljivk


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Uporabimo funkcijo unchop() za razširitev podatkov tabele bolezni_chop
# Chopped tabelo lahko razširimo za posamezno ali vse spremenljivke
head(unchop(bolezni_chop, c(Leto)), 3)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Razširitev glede na vse spremenljivke povrne tabelo  osnovno stanje
head(unchop(bolezni_chop, c(Leto, HepatitisB, HIV, Polio)), 3)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo pack() sesedemo stolpce HepatitisB, HIV in Polio v stolpec-tabelo "Bolezni". 
# Novonastalo tabelo poimenujemo bolezni_pack.
bolezni_pack <- pack(bolezni, Bolezen = c(HepatitisB, HIV, Polio))
# Preverimo tabelo
head(bolezni_pack, 3)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Ustvarili smo stolpec "Bolezen", ki ima eno vrstico za vsako državo in leto, 
# element vrstice pa je tabelo dimenzij 6x3
# Preverimo


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo colSums() izračunajmo vsoto 4. do 30. stolpca tabele zivDoba
colSums(zivDoba[, 4:13])


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# V stolpcih z manjkajočimi vrednostmi funkcija ne more izračunati vsote in vrne NA vrednost. 
# Obnašanje funkcije lahko prilagodimo z nastavitvijo parametra na.rm (na remove) na TRUE
colSums(zivDoba[, 4:13], na.rm = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo colMeans() izračunajmo vsoto 4. do 30. stolpca tabele zivDoba. 
# Funkcija naj odstani manjajoče vrednosti
colMeans(zivDoba[, 4:13], na.rm = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo quantile() izračunajmo kvantile z verjetnostjo 
# 0.05 in 0.95 spremenljivke Umrljivost v tabeli zibDoba
quantile(zivDoba$Umrljivost, probs = c(0.05, 0.95))


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo aggregate izračunajmo povprečno umrljivost v posamezni državi
aggregate(zivDoba$Umrljivost ~ zivDoba$Drzava, FUN="mean")

## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo aggregate izračunajmo še standardni odklon umrljivost glede na status države in leto
povpUmr_stLeto <- aggregate(zivDoba$Umrljivost ~ zivDoba$Status + zivDoba$Leto, FUN = "sd")
head(povpUmr_stLeto)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Najprej naložimo paket dplyr
library(dplyr)
# S funkcijo group_by() združimo podatke tabele zivDoba glede na leto in status 
# ter jih shranimo v zivDoba_group
zivDoba_group <- group_by(zivDoba, Leto, Status)
# Preverimo nastalo tabelo
head(zivDoba_group)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# S funkcijo summarize() izračunajmo povprećje in standardni odklon spremenljivk 
# Umrljivost in Alkohol v zivDoba_group. Shranimo podatke v zivDoba_sum
zivDoba_sum <- summarize(zivDoba_group, 
                         povpUmr = mean(Umrljivost), 
                         sdUmr = sd(Umrljivost), 
                         povpAlk = mean(Alkohol), 
                         sdAlk = sd(Alkohol))
head(zivDoba_sum)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# V enem koraku v zivDoba ustvarimo skupine glede na status in leto 
# ter izračunamo povprečje spremenljivke HepatitisB
zivDoba_hep <- zivDoba %>% group_by(Status, Leto) %>% 
                           summarize(povpHep = mean(HepatitisB, na.rm=T))
head(zivDoba_hep)


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# V enem koraku tabelo bolezni preuredimo s funkcijo pivot_longer(), 
# da bo imela vrednosti za vse bolezni v enem stolpcu (in imena bolezni v drugem), 
# združimo po državi in bolezni in izračunajmo povprečje
# Ponovno ustvarimo tabelo bolezni in naložimo tidyr, če jih nimamo več v delovnem okolju
bolezni_povp <- bolezni %>% pivot_longer(cols = c(HepatitisB, HIV, Polio), names_to = "Bolezen") %>% 
                            group_by(Drzava, Bolezen) %>% 
                            summarize(povpVrednost = mean(value))
head(bolezni_povp)

