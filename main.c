#include <stdio.h>      
#include <stdlib.h>   
#include <string.h> 
#include "kdtree.h"  
#include "fcm.h"    
#include "dbscan.h"   

int main(int argc, char *argv[]) {
    if (argc < 3) {
        return 1; 
    } //проверка

    int n_t, razm; // переменные для хранения количества точек и размерности пространства
    Point *pts = read_points(argv[1], &n_t, &razm); // Читаем точки из CSV файла
    if (!pts) {  printf("Error: cannot read file %s\n", argv[1]); return 1; } // ошибка
    printf("Loaded %d points of dimension %d\n", n_t, razm); // успешная загрузка


    KDNode *drv = NULL; // Указатель на корень дерева 
    for (int i = 0; i < n_t; i++) { drv = kd_insert(drv, pts[i], 0); } // вставляем все точки в дерево


    if (strcmp(argv[2], "-kd_insert") == 0) { // Проверка команды

        if (argc < 4) { //проверка аргумента
            printf("Error: specify coordinates\n");
            return 1; 
        } 

        double *crd = (double*)malloc(razm * sizeof(double)); // массив под координаты
        char *tok = strtok(argv[3], ",");   // Берем первое число из строки параметров
        for (int i = 0; i < razm && tok; i++) { crd[i] = atof(tok); tok = strtok(NULL, ","); } // Парсим все координаты
        Point nov = create_point(crd, razm); // Создаем структуру новой точки
        drv = kd_insert(drv, nov, 0);  // Вставляем её в построенное дерево
        printf("Point added\n"); // Уведомляем пользователя
        free(crd);                     
    }
    else if (strcmp(argv[2], "-kd_nearest") == 0) { // Проверка команды
        if (argc < 4) { 
            printf("Error: specify coordinates\n");
            return 1;
        } // Проверяем аргументы
        double *t_crd = (double*)malloc(razm * sizeof(double)); // память под целевую точку
        char *tok = strtok(argv[3], ",");   // Парсим строку с координатами цели
        for (int i = 0; i < razm && tok; i++) { t_crd[i] = atof(tok); tok = strtok(NULL, ","); } // Заполняем массив
        Point cel = create_point(t_crd, razm); // Создаем объект точки-цели
        double luch_d = 1e9; 
        KDNode *luch = NULL; // Инициализируем рекордное расстояние (очень большое) и указатель
        kd_find_nearest(drv, cel, 0, &luch, &luch_d); // Запускаем поиск в дереве
        if (luch) { // Если узел найден
            printf("Nearest point: ");
            for (int i = 0; i < razm; i++) { printf("%.2f ", luch->point.coords[i]); } // Выводим координаты найденной точки
            printf("\nDistance: %.2f\n", luch_d); // Выводим расстояние до неё
        }
        free_point(&cel); 
        free(t_crd);   
    }
    else if (strcmp(argv[2], "-cmeans") == 0) { // Проверка команды
        if (argc < 4) {
            printf("Error: specify number of clusters\n");
            return 1;
         } // Проверяем аргумент
        int k_t = atoi(argv[3]);// Преобразуем строку в число кластеров
        FCM *fc = fcm_create(pts, n_t, k_t, 2.0); // Создаем структуру FCM с параметром нечеткости 2.0
        fcm_run(fc, pts, 100, 0.001);  // Запускам алгоритм (макс 100 шагов, точность 0.001)
        fcm_print(fc, pts); 
        fcm_free(fc);                      
    }
    else if (strcmp(argv[2], "-dbscan") == 0) { // Проверка команды
        if (argc < 4) { 
            printf("Error: specify eps,minPts parameters\n");
            return 1; 
        } // Проверяем аргумент
        char *tok = strtok(argv[3], ","); double eps = atof(tok); // Парсим радиус окрестности
        tok = strtok(NULL, ","); int min_t = tok ? atoi(tok) : 3; // Парсим minPts 
        DBSCAN *db = dbscan_run(pts, n_t, eps, min_t); // Запускам алгоритм кластеризации
        dbscan_print(db, pts);
        dbscan_free(db);                   
    }
    else { 
        fprintf(stderr, "Unknown command: %s\n", argv[2]);// Если команда не распознана
    } 

    kd_free(drv);                
    for (int i = 0; i < n_t; i++) { 
        free_point(&pts[i]); 
    } 
    free(pts); return 0;               
}