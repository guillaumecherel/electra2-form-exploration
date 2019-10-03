# Explorer les formes

Electra 2 est une sculpture lumineuse cinétique. Un grand bras métallique de 2 mètres de long, fixé au mur en son centre, tourne entrainé par un moteur. Chaque moitié porte un autre bras motorisé, plus court. Tous les trois tiennent, à chaque extrêmité, une petite ampoule. On en joue, presque comme d'une marionnette, en actionnant du bout des doigts les commandes électriques qui contrôlent la vitesse des trois moteurs, leur sens de rotation et l'intensité des lumières. 

Dans l'obscurité, quand les ampoules s'allument et que les moteurs s'animent, des traces lumineuses s'impriment dans le noir. Ce sont d'abord des cercles, mais bien vite des boucles, des pointes, voire des lignes droites qui s'assemblent, se succèdent, se répètent et composent des rythmes et des formes surprenantes et captivantes.

La surprise a allumé notre curiosité et nous avons voulu chercher les formes que la sculpture pouvait dessiner. Les résultats sont exposés ici: une affiche qui donne une vue d'ensemble statique des formes que nous avons découvertes et une interface interactive pour les voir s'animer. La suite explique notre démarche.

Il y a quelques années à l'Institut des Systèmes Complexes, nous avons développé un algorithme pour rechercher les différents résultats que peut donner un programme informatique dont l'exécution dépend de données d'entrée. L'objectif était de rechercher ses comportements inattendus, par exemple pour repérer des erreurs. Nous avons eu l'idée de l'utiliser pour trouver des formes dessinées par Electra 2.

Pour employer cet algorithme, nous devions donc modéliser Electra 2 sous forme d'un programme informatique. Il calcule les trajectoires des ampoules par de simples formules géométriques, à partir des longueurs des bras et de leurs vitesses de rotation.

Pour pouvoir distinguer les formes rares des plus fréquentes, l'algorithme a besoin d'une description succinte, en quelques nombres, de chaque forme dessinée. La plus grande difficulté a été de trouver une telle description qui permettait de distinguer les formes qui nous semblent différentes à l'œil nu.

Nos premières observations de la sculpture nous avait donné une première idée des caractéristiques que l'on pouvait attendre: boucles, segments rectilignes, pointes, répétitions. Nous avons choisi de décrire les formes par leur courbure, le nombre de points singuliers (les points de la trajectoire où la vitesse d'une lumière est nulle, comme lors de l'apparition de pointes) et le temps que mettent les trajectoires à se répéter. 

Avec ces mesures, notre algorithme nous a permis de rechercher les formes qui combinaient ces différentes caractéristiques. Avec l'interface interactive, vous pouvez explorer les formes en spécifiant des valeurs qui les décrivent. Pour créer l'affiche, nous avons regroupé les formes similaires. Chaque trajectoire représente un groupe.

