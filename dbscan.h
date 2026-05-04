#ifndef DBSCAN_H
#define DBSCAN_H

#include "kdtree.h"

// Статусы точек (до кластеризации)
#define UNCLASSIFIED -1 // еще не обработана
#define NOISE 0// шум (не принадлежит ни одному кластеру)

// Структура для результатов DBSCAN
typedef struct {
    int *labels; // метка кластера для каждой точки
    int n_points;// сколько точек
    int n_clusters;// сколько кластеров нашли
} DBSCAN;


// главный цикл
DBSCAN* dbscan_run(Point *points, int n_points, double eps, int min_pts);
// Печать результатов
void dbscan_print(DBSCAN *db, Point *points); 
//очистка
void dbscan_free(DBSCAN *db);   

#endif