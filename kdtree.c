#include "kdtree.h"   
#include <stdio.h>         
#include <stdlib.h>   
#include <math.h>   
#include <string.h>  


// Функция создания новой точки
Point create_point(double *crd, int razm) {
    Point tch;                         
    tch.razmer = razm;// размерность передаем( 3д, 2д )
    tch.coords = (double*)malloc(razm * sizeof(double)); //массив под координаты
    for (int i = 0; i < razm; i++) {//проход по всем осям пространства( x- первое число, y - второе и тд)
        tch.coords[i] = crd[i];//Копируем значение координаты в массив
    }
    return tch;
}

// Функция очистки памяти точки
void free_point(Point *tch) { 
    if (tch->coords) {// проверка выделения памяти
        free(tch->coords);//само освобождение
        tch->coords = NULL;
    }
}
// Функция создания узла k-d дерева
KDNode* kd_create_node(Point tch, int os) { 
    KDNode *uzel = (KDNode*)malloc(sizeof(KDNode)); //память под узел
    uzel->point.razmer = tch.razmer;// размерность копируем в узел
    uzel->point.coords = (double*)malloc(tch.razmer * sizeof(double)); // выделяем НОВУЮ память для копии
    memcpy(uzel->point.coords, tch.coords, tch.razmer * sizeof(double)); // копируем данные в новую память - то есть уже прям в узел
    uzel->levo = NULL;        
    uzel->pravo = NULL;                 
    uzel->axis = os;                   
    return uzel;                           
}

// Рекурсивная вставка точки
KDNode* kd_insert(KDNode *koren, Point tch, int glub) { //корень - узел текущ, tch - наша точкаm 
    if (koren == NULL) { // Если дошли до пустого места 
        return kd_create_node(tch, glub % tch.razmer); // создаем новый узел
    }
    int os = koren->axis; // берем ось
    if (tch.coords[os] < koren->point.coords[os]) { // если координата новой точки меньше
        koren->levo = kd_insert(koren->levo, tch, glub + 1); //идем в левое поддерево
    } else { //если координата больше или равна
        koren->pravo = kd_insert(koren->pravo, tch, glub + 1); //идем в правое поддерево
    }
    return koren;                       
}

// Расчет расстояния между двумя точками по Евклиду
double euclidean_distance(Point t1, Point t2) { 
    if (t1.razmer != t2.razmer) {// проверка на размерность
        return -1;
     } 
    double r = 0.0;//для суммы
    for (int i = 0; i < t1.razmer; i++) {//проход по всем координатам точек
        double raz = t1.coords[i] - t2.coords[i]; // считаем разницу по текущ оси
        r += raz * raz; //добавляем квадрат разницы к общей сумме
    }
    return sqrt(r); //евклидово расстояние
}

// Поиск ближайшего соседа
KDNode* kd_find_nearest(KDNode *koren, Point cel, int glub, KDNode **luch, double *luch_d) { 
    if (koren == NULL) { //если узел пустой, возвращаем лучший из найденных ранее
        return *luch;     
    } 

    double r = euclidean_distance(koren->point, cel); // считаем расстояние
    if (r < *luch_d) { //Если нашли точку ближе текущ рекорда
        *luch_d = r;//обновляем переменную с лучшим расстоянием
        *luch = koren;//запоминаем указатель на этот узел как лучший
    }
    int os = koren->axis;//	Ось, по которой разделено пространство в текущем узле
    KDNode *vper, *vtor;//Указатели на "перспективную" и "вторую" ветки

    if (cel.coords[os] < koren->point.coords[os]) {//цель левее -> сначала идем влево
        vper = koren->levo; 
        vtor = koren->pravo; 
    } 
    else { //цель правее -> сначала идем вправо
        vper = koren->pravo; 
        vtor = koren->levo; 
    }
    kd_find_nearest(vper, cel, glub + 1, luch, luch_d); //Рекурсивно ищем в первой (ближней) ветке
    double os_r = fabs(cel.coords[os] - koren->point.coords[os]); //Считаем расстояние до разделяющей плоскости
    if (os_r < *luch_d) {  // Если плоскость ближе рекорда, там может быть ближайшая точка
        kd_find_nearest(vtor, cel, glub + 1, luch, luch_d); // Проверяем и вторую ветку
    }
    return *luch;           
}

// Рекурсивное удаление дерева из памяти
void kd_free(KDNode *koren) { 
    if (koren == NULL) return; 
    kd_free(koren->levo);
    kd_free(koren->pravo);          
    free_point(&(koren->point)); 
    free(koren); 
}


// Чтение точек из CSV файла
Point* read_points(char *fname, int *n_t, int *razm) { 
    FILE *f = fopen(fname, "r");      
    if (!f) {
        printf("Не могу открыть файл %s\n", fname); 
        return NULL; 
    } 
    char str[1024]; 
    *n_t = 0; // буфер для строки и счетчик колва точек
    while (fgets(str, sizeof(str), f)) (*n_t)++; // читаем файл построчно, считая строки

    rewind(f);//возвращаем курсор чтения в начало файла

    fgets(str, sizeof(str), f);
     *razm = 1; // читаем первую строку, начальная размерность 1

    for (char *c = str; *c; c++) if (*c == ',') (*razm)++; // считаем запятые, чтобы понять размерность

    rewind(f);// снова возвращаем курсор в начало

    Point *pts = (Point*)malloc((*n_t) * sizeof(Point)); // массив структур Point
    for (int i = 0; i < *n_t; i++) {
        fgets(str, sizeof(str), f);//читаем текущую строку с координатами
        pts[i].razmer = *razm;// Записываем размерность в структуру точки
        pts[i].coords = (double*)malloc((*razm) * sizeof(double)); // память под массив координат
        char *tok = strtok(str, ",\n");// берем первое число до запятой или конца строки

        for (int j = 0; j < *razm && tok; j++) { 
            pts[i].coords[j] = atof(tok);// Преобразуем строку в double и записываем
            tok = strtok(NULL, ",\n");//Берем следующее число
        }
    }
    fclose(f); 
    return pts; 
}
