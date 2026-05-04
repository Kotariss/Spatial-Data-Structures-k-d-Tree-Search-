#include "dbscan.h" 
#include <stdio.h>  
#include <stdlib.h>       
#include <string.h>       
#include <math.h>        

// Поиск соседей в радиусе eps
int* find_neighbors(Point *pts, int n_t, Point p, double eps, int *n_sos) {//n_t - колво точек

    int *sos = (int*)malloc(n_t * sizeof(int)); //массив максимальной длины 
    *n_sos = 0; //сначало соседей 0

    for (int i = 0; i < n_t; i++) {
        if (euclidean_distance(pts[i], p) <= eps) { //если расстояние до точки p меньше радиуса
            sos[(*n_sos)++] = i; // Записываем индекс точки и увеличиваем счетчик
        }
    }
    return sos; 
}

 // Расширение кластера
void expand_cluster(Point *pts, int *met, int n_t, int idx, int id_kl, double eps, int min_t) {
    int n_sos; int *sos = find_neighbors(pts, n_t, pts[idx], eps, &n_sos); //ищем соседей начальной точки
    if (n_sos < min_t) {//Если соседей мало, помечаем как шум и выходим
        met[idx] = NOISE;
        free(sos); 
        return; 
    } 
    for (int i = 0; i < n_sos; i++) {// Иначе помечаем всех найденных соседей меткой кластера
         met[sos[i]] = id_kl; 
        } 
    int *sem = (int*)malloc(n_sos * sizeof(int)); // Создаем очередь для точек, которые нужно проверить
    memcpy(sem, sos, n_sos * sizeof(int));// Копируем индексы соседей в очередь
    int n_sem = n_sos; // текущ размер очереди
    for (int i = 0; i < n_sem; i++) { 
        if (sem[i] == idx) { 
            sem[i] = sem[--n_sem]; 
            break;
        } // Заменяем на последний элемент и уменьшаем размер
    }
    while (n_sem > 0) { // Пока очередь точек для проверки не пуста

        int tek = sem[--n_sem]; // Берем последнюю точку из очереди(потому что используем стек (LIFO), но для дбскана не важно, просто нам так легче)*
        int n_tek_sos; int *tek_sos = find_neighbors(pts, n_t, pts[tek], eps, &n_tek_sos); // Ищем соседей этой точки

        if (n_tek_sos >= min_t) { // Если это тоже "ядерная" точка (ну или центральная)
            for (int i = 0; i < n_tek_sos; i++) { //проходим по всем её соседям
                int nb = tek_sos[i];// Индекс соседа
                if (met[nb] == UNCLASSIFIED || met[nb] == NOISE) { // Если сосед еще не в кластере или помечен шумом
                    if (met[nb] == UNCLASSIFIED) { // Если совсем не обработан
                        
                        sem = (int*)realloc(sem, (n_sem + 1) * sizeof(int)); // Увеличиваем память очереди
                        sem[n_sem++] = nb; // Добавляем индекс соседа в очередь для дальнейшего расширения
                        //тем не менее мы не будем его проверять тут, тк  n_sem > 0 - он банально уйдет в новый кластер*
                    }
                    met[nb] = id_kl;// В любом случае присваиваем ему текущую метку кластера
                }
            }
        }
        free(tek_sos); 
    }
    free(sos);
    free(sem);  
}


//DBSCAN естественно
DBSCAN* dbscan_run(Point *pts, int n_t, double eps, int min_t) { 
    DBSCAN *db = (DBSCAN*)malloc(sizeof(DBSCAN)); // результаты

    db->n_points = n_t; 
    db->labels = (int*)malloc(n_t * sizeof(int)); //массив меток для каждой точки

    for (int i = 0; i < n_t; i++) { // Изначально все точки не классифицированы
        db->labels[i] = UNCLASSIFIED; 
    } 
    int id_kl = 0;// Счетчик найденных кластеров
    for (int i = 0; i < n_t; i++) { // Проходим по всем точкам датасета
        if (db->labels[i] == UNCLASSIFIED) { // Если точка еще не тронута алгоритмом
            expand_cluster(pts, db->labels, n_t, i, ++id_kl, eps, min_t); // Запускаем расширение нового кластера
        }
    }
    db->n_clusters = id_kl; //итоговое число
    return db;
}

// Вывод 
void dbscan_print(DBSCAN *db, Point *pts) { 
    printf("\nDBSCAN RESULTS\nFound clusters: %d\n", db->n_clusters); 

    int *razm_k = (int*)calloc(db->n_clusters + 1, sizeof(int));
    int shum = 0;

    for (int i = 0; i < db->n_points; i++) { 
        if (db->labels[i] == NOISE) { // Если метка NOISE, увеличиваем счетчик шума
            shum++; 
        } 
        else {// Иначе увеличиваем счетчик точек в соответствующем кластере
            razm_k[db->labels[i]]++; 
        } 
    }
    printf("Noise points: %d\nCluster sizes:\n", shum); // Выводим количество выбросов
    for (int i = 1; i <= db->n_clusters; i++) {
        printf("Cluster %d: %d points\n", i, razm_k[i]);
    }
    printf("\nPoint assignment:\n");
    for (int i = 0; i < db->n_points; i++) {

        printf("Point %d (", i);// Выводим номер точки

        for (int d = 0; d < pts[i].razmer; d++) { // Выводим её координаты

            printf("%.2f ", pts[i].coords[d]);
        }
        printf( ") -> cluster %d\n", db->labels[i]); // выводим итоговую метку
    }
    free(razm_k);
}

void dbscan_free(DBSCAN *db) {  
    free(db->labels);   
    free(db);       
}
