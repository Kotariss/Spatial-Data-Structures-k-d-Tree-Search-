#ifndef KDTREE_H 
#define KDTREE_H

// Структура для хранения точки
typedef struct {
    double *coords;// координат x,y,z
    int razmer;// размерность пространства 3D-3,2D-2
} Point;

// Структура для узла k-d дерева
typedef struct KDNode {
    Point point;// точка, хранящаяся в этом узле
    struct KDNode *levo;// указатель на левое поддерево
    struct KDNode *pravo;// указатель на правое поддерево
    int axis; // по какой оси делим (0 для x, 1 для y, 2 для z)
} KDNode;

//все функции
Point* read_points(char *filename, int *n, int *d); // читаем точки из файла

Point create_point(double *coords, int razmer);// создаем новую точку
KDNode* kd_create_node(Point p, int axis);// создаем узел дерева

KDNode* kd_insert(KDNode *root, Point p, int depth);// вставляем точку в дерево

double euclidean_distance(Point a, Point b);// считаем расстояние между точками

KDNode* kd_find_nearest(KDNode *root, Point target,int depth, KDNode **best, double *best_dist);// ищем ближайшую точку

void free_point(Point *p);// удаляем точку из памяти
void kd_free(KDNode *root); // удаляем все дерево

#endif