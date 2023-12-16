#include <stdio.h>
#include "vertices.cpp"


void fill_vertices(Vertices &ver, int n)
{
    double x[3] = {0.1, 0.2, 0.3};
    for (int i=0;i<n;i++) {
        ver.append_node(x, 100+i, i);
        x[0]++; x[1]++; x[2]++;
    }
}

void print_vertices(Vertices &ver)
{
    printf("%d vertices:\n", ver.size);
    for (int i=0;i<ver.size;i++) {
        double *x = ver.pos + 3*i;
        printf("id=%2d, x=%4.1f, y=%4.1f, z=%4.1f, idtri=%3d, sign=%2d\n",
               i, x[0], x[1], x[2], ver.idtri[i], ver.sign[i]);
    }
}

void test_renumber1(Vertices &ver)
{
    int idx[3] = {5, 6, 0};
    printf("new indices: %d %d %d\n", idx[0], idx[1], idx[2]);
    int invalid = ver.renumber_indices(3, idx);
    printf("new indices: %d %d %d\n", idx[0], idx[1], idx[2]);
    printf("#invalid = %d\n", invalid);
}

void test_renumber2(Vertices &ver)
{
    int idx[3] = {2, 3, 4};
    printf("new indices: %d %d %d\n", idx[0], idx[1], idx[2]);
    int invalid = ver.renumber_indices(3, idx);
    printf("new indices: %d %d %d\n", idx[0], idx[1], idx[2]);
    printf("#invalid = %d\n", invalid);
}

int main()
{
    Vertices ver;
    ver.allocate(10);
    fill_vertices(ver, 9);
    print_vertices(ver);
    int flag[9] = {1, 0, 0, 1, 0, 1, 1, 0, 0};
    int cnt = ver.remove_if_zero(flag);
    printf("remove %d nodes\n", cnt);
    print_vertices(ver);

    test_renumber1(ver);
    test_renumber2(ver);


}