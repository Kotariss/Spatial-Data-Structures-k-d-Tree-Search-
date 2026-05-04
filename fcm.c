#include "fcm.h"      
#include <stdio.h>       
#include <stdlib.h>     
#include <math.h>      
#include <time.h>     

// Создание структуры FCM
// pts — точки данных, n_t — их число, n_k — кластеров, fuz — нечёткость
FCM* fcm_create(Point *pts, int n_t, int n_k, double fuz) {
    FCM *fc = (FCM*)malloc(sizeof(FCM));      
    fc->n_points = n_t; 
    fc->n_clusters = n_k; 
    fc->fuzziness = fuz;

    fc->centroids = (Point*)malloc(n_k * sizeof(Point));// Массив центров
    for (int i = 0; i < n_k; i++) { // для каждого кластера
        // память под координаты центра
        fc->centroids[i].coords = (double*)malloc(pts[0].razmer * sizeof(double));
        fc->centroids[i].razmer = pts[0].razmer; // размерность
    }

    // Выделяем память под матрицу принадлежности [точки × кластеры]
    fc->membership = (double**)malloc(n_t * sizeof(double*));
    for (int i = 0; i < n_t; i++) { // Для каждой точки
        fc->membership[i] = (double*)malloc(n_k * sizeof(double)); // строка матрицы
        
        double row_sum = 0.0; // Сумма для нормировки
        for (int j = 0; j < n_k; j++) { // Для каждого кластера
            // Случайное число [0.5, 1.5] чтобы избежать нулей
            fc->membership[i][j] = 0.5 + (double)rand() / RAND_MAX;
            row_sum += fc->membership[i][j];// Накопление суммы
        }
        for (int j = 0; j < n_k; j++) {// Нормировка строки
            fc->membership[i][j] /= row_sum; // делим на сумму 
        }
    }
    return fc;
}

// Обновление матрицы принадлежности
void fcm_update_membership(FCM *fc, Point *pts) { 
    double m = fc->fuzziness;// нечеткость
    double st = 2.0 / (m - 1.0);//показатель степени для формулы

    for (int i = 0; i < fc->n_points; i++) {
        double summ = 0.0; //накопитель знаменателя

        double *rass = (double*)malloc(fc->n_clusters * sizeof(double));// Массив расстояний до центров
        int zero_cluster = -1;// Флаг

        for (int j = 0; j < fc->n_clusters; j++) { 
            // Считаем евклидово расстояние
            rass[j] = euclidean_distance(pts[i], fc->centroids[j]);
            if (rass[j] == 0.0) { // Если точка == центр
                zero_cluster = j;// запоминаем индекс
                break;
            }
        }

        if (zero_cluster != -1) {// если нашли точное совпадение
            for (int j = 0; j < fc->n_clusters; j++) { 
                // 1.0 своему, 0.0 чужим
                fc->membership[i][j] = (j == zero_cluster) ? 1.0 : 0.0;
            }
            free(rass);                         
            continue;   
        }

        for (int j = 0; j < fc->n_clusters; j++) { // Считаем знаменатель
            // Суммируем (1/d)^st для нормировки
            summ += pow(1.0 / rass[j], st);
        }
        
        for (int j = 0; j < fc->n_clusters; j++) { // Вычисляем доли
            // Формула принадлежности: числитель / знаменатель
            fc->membership[i][j] = pow(1.0 / rass[j], st) / summ;
        }
        free(rass);                      
    }
}

// Обновление центров кластеров 
void fcm_update_centers(FCM *fc, Point *pts) { 
    double m = fc->fuzziness;//нечеткость
    int razm = pts[0].razmer; //размерность пространства

    for (int j = 0; j < fc->n_clusters; j++) {  
        //по формуле
        // числитель
        double *chisl = (double*)calloc(razm, sizeof(double));
        double znam = 0.0; //знаменатель

        for (int i = 0; i < fc->n_points; i++) { 
            // Вес точки: u^m
            double u = pow(fc->membership[i][j], m);
            for (int d = 0; d < razm; d++) { 
                // Накопление: вес × координата точки
                chisl[d] += u * pts[i].coords[d];
            }
            znam += u;//накопление суммы весов
        }
        
        for (int d = 0; d < razm; d++) {
            // Формула: числитель / знаменатель (с защитой от деления на 0)
            fc->centroids[j].coords[d] = (znam > 1e-10) ? chisl[d] / znam : 0.0;
        }
        free(chisl);
    }
}



// САМ ЦИКЛ АЛГОРИТМА!!!
void fcm_run(FCM *fc, Point *pts, int max_it, double epsil) {

    fcm_update_centers(fc, pts);// вычисляем центры из случайной матрицы U

    double pred_osh = 1e9;// ошибка

    for (int iter = 0; iter < max_it; iter++) {

        Point *star = (Point*)malloc(fc->n_clusters * sizeof(Point)); //для сохранения старых центров
        for (int j = 0; j < fc->n_clusters; j++) {
            // Глубокое копирование точки
            star[j] = create_point(fc->centroids[j].coords, pts[0].razmer);
        }

        fcm_update_membership(fc, pts);//обновляем матрицу U
        fcm_update_centers(fc, pts); //обновляем центры V

        double osh = 0.0; //cуммарный сдвиг центров
        for (int j = 0; j < fc->n_clusters; j++) { 
            // Считаем расстояние сдвига центра
            osh += euclidean_distance(fc->centroids[j], star[j]);
            free_point(&star[j]);
        }
        free(star); 

        // Критерий остановки: центры почти не двигаются
        if (osh < epsil) {
            // Вывод сообщения о сходимости
            printf("Converged at iteration %d (shift: %e < %e)\n", iter, osh, epsil);
            break;                         
        }
        pred_osh = osh;                    
    }
}

// Вывод результатов
void fcm_print(FCM *fc, Point *pts) { 
    printf("\nFCM RESULTS\n");    
    // Параметры алгоритма
    printf("Clusters: %d, Points: %d, Fuzziness: %.2f\n", 
           fc->n_clusters, fc->n_points, fc->fuzziness);
    
    printf("Cluster centers:\n");           
    for (int j = 0; j < fc->n_clusters; j++) { 
        printf("  K%d: [", j);// Номер кластера
        for (int d = 0; d < pts[0].razmer; d++) {
            // Печать координаты с запятой или без
            printf("%.3f%s", fc->centroids[j].coords[d], (d < pts[0].razmer-1) ? ", " : "");
        }
        printf("]\n");                          
    }
    
    printf("\nPoint membership (first 10):\n");
    // Ограничиваем вывод 10 точками
    int limit = (fc->n_points < 10) ? fc->n_points : 10;
    for (int i = 0; i < limit; i++) { 
        printf("  P%d: [", i);// Номер точки
        for (int j = 0; j < fc->n_clusters; j++) { 
            // Печать степени принадлежности
            printf("%.3f%s", fc->membership[i][j], (j < fc->n_clusters-1) ? ", " : "");
        }
        printf("]\n");                        
    }
    // Если точек больше 10 — сообщаем
    if (fc->n_points > 10) printf("  ... and %d more points\n", fc->n_points - 10);
}

// Очистка памяти
void fcm_free(FCM *fc) { 
    for (int j = 0; j < fc->n_clusters; j++) {  
        free_point(&fc->centroids[j]);          
    }
    free(fc->centroids);                        
    
    for (int i = 0; i < fc->n_points; i++) {    
        free(fc->membership[i]);               
    }
    free(fc->membership);                     
    free(fc);                                
}