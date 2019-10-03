# Explorer les formes

Electra 2 est une sculpture lumineuse cinétique. Un grand bras métallique de 2 mètres de long, fixé au mur en son centre, tourne entrainé par un moteur. Chaque moitié porte un autre bras motorisé, plus court. Tous les trois tiennent, à chaque extrêmité, une petite ampoule. On en joue, presque comme d'une marionnette, en actionnant du bout des doigts les commandes électriques qui contrôlent la vitesse des trois moteurs, leur sens de rotation et l'intensité des lumières. 

Dans l'obscurité, quand les ampoules s'allument et que les moteurs s'animent, des traces lumineuses s'impriment dans le noir. Ce sont d'abord des cercles, mais bien vite des boucles, des pointes, voire des lignes droites qui s'assemblent, se succèdent, se répètent et composent des rythmes et des formes surprenantes et captivantes.

La surprise a allumé notre curiosité et nous avons voulu chercher les formes que la sculpture pouvait dessiner. Les résultats sont exposés ici: une affiche qui donne une vue d'ensemble statique des formes que nous avons découvertes et une interface interactive pour les voir s'animer. La suite explique notre démarche.

Il y a quelques années à l'Institut des Systèmes Complexes, nous avons développé un algorithme pour rechercher les différents résultats que peut donner un programme informatique dont l'exécution dépend de paramètres en entrée. L'objectif était de rechercher ses comportements inattendus, par exemple pour repérer des erreurs. Nous avons eu l'idée de l'utiliser pour trouver des formes dessinées par Electra 2.

Pour employer cet algorithme, nous devions donc modéliser Electra 2 sous forme d'un programme informatique. Il calcule les trajectoires des ampoules par de simples formules géométriques, à partir des longueurs des bras et de leurs vitesses de rotation.

Pour pouvoir distinguer les formes rares des plus fréquentes, l'algorithme a besoin d'une description succinte, en quelques nombres, de chaque forme dessinée. La plus grande difficulté a été de trouver une telle description qui permettait de distinguer les formes qui nous semblent différentes à l'œil nu.

Nos premières observations de la sculpture nous avait donné une première idée des caractéristiques que l'on pouvait attendre: boucles, segments rectilignes, pointes, répétitions. Nous avons choisi de décrire les formes en mesurant pour chaque trajectoire (une trajectoire correspondant à une ampoule) la courbure, le nombre de ralentissements (les points de la trajectoire où la vitesse d'une ampoule passe sous un seuil fixé) et la longueur du motif (c'est-à-dire le temps que met une trajectoire à se répéter). Pour simplifier le travail, nous avons mesuré ces caractéristiques pour seulement deux trajectoires parmi les six possibles, une sur chacun des deux petits bras.

Avec ces mesures, notre algorithme nous a permis de rechercher les formes qui combinaient ces différentes caractéristiques. L'affiche représente un échantillon de celles que nous avons découvertes. Chaque dessin représente les deux trajectoires que nous avons étudiées et elles sont aussi longues que nécessaires pour que les motifs se répètent au moins 6 fois. La transparence reflète la vitesse des lumières.

Avec l'interface interactive, vous pouvez voir ces formes s'animer. Essayez de retrouver celles de l'affiche en contrôlant les caractéristiques de chaque trajectoire: courbure, ralentissements et longueur de motif. 



# Directory structure

- formExploration/
  - Openmole/electra: contain .oms files to run the model in openmole and perform different analyses/explorations : directSampling, NSGA (is it possible to produce a given trajectory with electra?), and PSE (explore the diversity in trajectories)  
  - R electra: We have also implemented the model in R.
    - R model: R code to plot trajectories (we indicate some parameters that produce simple figures) and to save them (to create a movie)
    - postProc_OM: R code to analyse results from openmole, sorted by explicit name directory (NSGA,PSE,direcSampling,1run). The final figure is obtained in the PSE folder.
  - scala: contain the scala code of the model. In particular there is a function to produce the stationary trajectory (also a function for the non stationary bu not used yet),  It also contains functions that compute measures on the trajectories (loop points, Moran index, density index, curvature and singular points.

