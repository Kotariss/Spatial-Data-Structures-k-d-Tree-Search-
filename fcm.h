#ifndef FCM_H
#define FCM_H

#include "kdtree.h"

// Структура для алгоритма fcm
typedef struct {
    Point *centroids; // центры кластеров
    double **membership;// матрица принадлежности [точка][кластер]
    int n_points;// сколько всего точек
    int n_clusters;// сколько кластеров
    double fuzziness;// "нечеткость" (обычно 2)
} FCM;


// Создаёт и инициализирует структуру FCM
FCM* fcm_create(Point *points, int n_points, int n_clusters, double fuzziness);

// Обновляет матрицу принадлежности U по текущим центрам
//uᵢⱼ = 1 / Σ(dᵢⱼ/dₖⱼ)^(2/(m-1)) формцла которую используем
void fcm_update_membership(FCM *fcm, Point *points);

// Обновляет центры кластеров V по матрице принадлежности
//vⱼ = Σ(uᵢⱼ^m * xᵢ) / Σ(uᵢⱼ^m) формцла которую используем
void fcm_update_centers(FCM *fcm, Point *points);

// главный цикл
void fcm_run(FCM *fcm, Point *points, int max_iter, double epsilon);

// Вывод результатов
void fcm_print(FCM *fcm, Point *points);

// Освобождает всю память
void fcm_free(FCM *fcm);

#endif