import numpy as np
import os
import scipy.spatial
from scipy.spatial import distance


def lee_coordenadas_atomo(linea):
    """Interpreta las coordenadas de una línea de un fichero PDB que empiece por
    ATOM (es un átomo)"""
    if linea.startswith('ATOM  '):
        x = float(linea[30:38])
        y = float(linea[38:46])
        z = float(linea[46:54])
        return [x, y, z]


def obtiene_coordenadas(estructura):
    """Lee una estructura (fichero PDB) y obtiene las coordenadas de los
    átomos que contiene"""
    coordenadas = []
    with open(estructura) as input:
        lineas = [linea.rstrip(os.linesep) for linea in input.readlines()]
        for linea in lineas:
            atomo = lee_coordenadas_atomo(linea)
            if atomo:
                coordenadas.append(atomo)
    return np.array(coordenadas)

# Descargando datos
!wget -nc https://gitlab.uoclabs.uoc.es/prog_bioinf/data/-/raw/master/Unidad3/1PPE_rec.pdb --no-check-certificate
!wget -nc https://gitlab.uoclabs.uoc.es/prog_bioinf/data/-/raw/master/Unidad3/1PPE_lig.pdb --no-check-certificate

# Coordenadas_1 tiene las coordenadas de la proteína A en formato numpy
coordenadas_1 = obtiene_coordenadas('1PPE_rec.pdb')

# Coordenadas_2 tiene las coordenadas de la proteína B en formato numpy
coordenadas_2 = obtiene_coordenadas('1PPE_lig.pdb')

# Código a completar:

# Pista: podéis utilizar la función scipy.spatial.distance.cdist y la función numpy.where.


def smaller_than(x,threshold):
  # np.where nos da las coordenadas del array de los valores del array en los que se cumple la condición.
  # voy a crear un array que contenga los índices de aquellos valores del array cumplan la condición (valor en x mayor que el threshold marcado).
  # este array se compone de tuplas, las cuales nos dan la ubicación de aquellos valores que reportan esa condición.
  # la fila se corresponde con la línea del documento en la que se cumple esa condición, cada línea se corresponde con un átomo de la proteína.
  index = np.where(x <= threshold)
  list_index = index[0].tolist() # convierto en lista la primera fila del array anterior que se corresponde con las filas(átomos) de la proteína que uso como referencia que cumplen la condición.
  return list_index # acabamos teniendo una lista con una serie de números que nos dice en qué filas del array se ha cumplido la condición y si aparece más de una vez significa que ese átomo ha cumplido varias veces esa condición.


def atomos_contacto(lista):
  # hacemos una función que itere sobre los elementos de una lista y los convierta en claves de un diccionario. Los elementos de esa lista son índices.
  # las veces que se repite ese elemento/clave en la lista serán los valores de esas claves en el diccionario.
  atomos_in_contact = {} # creamos diccionario
  for indice in lista: # iteramos sobre los elementos de la lista
    if indice not in atomos_in_contact.keys(): # si ese elemento no está entre las claves del diccionario:
      atomos_in_contact[indice] = 1 # lo agregamos como clave y le damos el valor de 1
    else: # si ese elemento de la lista ya está como clave en el diccionario:
      atomos_in_contact[indice] += 1 # sumamos 1 al valor de esa clave cada vez que eso ocurra.
  return len(atomos_in_contact.keys()) # la función nos devuelve la longitud del número de claves en el diccionario 'atomos_in_contact'.

mod_atomos_A_B = distance.cdist(coordenadas_1, coordenadas_2, 'euclidean') # obtenemos el array de las distancias existentes entre todos los átomos de A y los de B, filas -> átomos de A, columnas -> átomos de B.


atomos_A_cercanos = smaller_than(mod_atomos_A_B,3) # obtenemos la lista de los índices de las filas (átomos de A) en las que se encuentran los valores que han cumplido
                                                   # la condición de que las distancias entre átomos de A y B son iguales o inferiores a 3Å.

atomos_A_B = atomos_contacto(atomos_A_cercanos) # creamos diccionario a partir de la lista creada con la función 'smaller_than'. Las claves son los índices de las filas (átomos) de la proteína A y valores las veces que se repiten.
                                                # la función 'atomos_contacto' nos devuelve el número de claves del diccionario, que serán el número de átomos de A en contacto con los de B.


mod_atomos_B_A = distance.cdist(coordenadas_2, coordenadas_1, 'euclidean') # obtenemos el array de las distancias existentes entre todos los átomos de B y los de A. filas -> átomos de B, columnas -> átomos de A.
                                                                           # las medidas son iguales a las obtenidas en 'mod_atomos_A_B', pero se invierten filas y columnas
                                                                           # porque cambia la proteína de referencia A<->B.

atomos_B_cercanos = smaller_than(mod_atomos_B_A,3) # obtenemos la lista de los índices de las filas (átomos de B) en las que se encuentran los valores que han cumplido
                                                   # la condición de que las distancias entre átomos de B y A son iguales o inferiores a 3Å.

atomos_B_A = atomos_contacto(atomos_B_cercanos) # creamos diccionario a partir de la lista creada con la función 'smaller_than'. Las claves son los índices de las filas (átomos) de la proteína B y valores las veces que se repiten.n.
                                                # la función 'atomos_contacto' nos devuelve el número de claves del diccionario, que serán el número de átomos de B en contacto con los de A.


print("Número de átomos de A en contacto con B: ", atomos_A_B)

print("Número de átomos de B en contacto con A: ", atomos_B_A)

# Finalmente, el número de átomos total será la suma de ambos:
print("Número total de átomos en contacto: ", atomos_A_B + atomos_B_A)
