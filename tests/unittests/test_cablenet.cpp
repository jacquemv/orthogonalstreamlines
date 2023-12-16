#include <stdio.h>
#include "cablenetwork.cpp"

#define VERBOSE 0

void print_cablenet(CableNetwork &cnet)
{
    if (!VERBOSE) return;
    printf("%d cables:\n", cnet.nc);
    for (int c=0;c<cnet.nc;c++) {
        printf("cable %d: ", c);
        for (int i=cnet.sep[c];i<cnet.sep[c+1];i++)
            printf("%d ", cnet.idx[i]);
        if (cnet.group_sep == c+1)
            printf(" /");
        printf("\n");
    }
}

void test_creation()
{
    CableNetwork cnet;
    cnet.allocate(5, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(3); cnet.append(1);
    cnet.new_cable();
    cnet.append(2); cnet.append(4); cnet.append(5);
    cnet.new_group();
    cnet.new_cable();
    cnet.append(9);
    print_cablenet(cnet);
    assert(cnet.size() == 6);
    assert(cnet.nc == 3);
    assert(cnet.group_sep == 2);
    assert((cnet.idx[0] == 3) && (cnet.idx[1] == 1) && (cnet.idx[2] == 2)
        && (cnet.idx[3] == 4) && (cnet.idx[4] == 5) && (cnet.idx[5] == 9));
    assert((cnet.sep[0] == 0) && (cnet.sep[1] == 2) && (cnet.sep[2] == 5) 
        && (cnet.sep[3] == 6));
}

void test_squeeze0()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    print_cablenet(cnet);
    int cnt = cnet.squeeze();
    print_cablenet(cnet);
    if (VERBOSE) printf("removed %d cables\n", cnt);
    assert(cnet.size() == 0);
    assert(cnet.nc == 0);
    assert(cnet.group_sep == 0);
    assert(cnet.sep[0] == 0);
}

void test_squeeze1()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.new_cable();
    cnet.append(3); cnet.append(1);
    cnet.new_cable();
    cnet.append(7);
    cnet.new_cable();
    cnet.new_group();
    cnet.new_cable();
    cnet.new_cable();
    cnet.append(2); cnet.append(4); cnet.append(5);
    cnet.new_cable();
    cnet.append(9);
    cnet.new_cable();
    print_cablenet(cnet);
    assert(cnet.size() == 7);
    assert(cnet.largest_index() == 9);
    int cnt = cnet.squeeze();
    print_cablenet(cnet);
    if (VERBOSE) printf("removed %d cables\n", cnt);
    assert(cnet.nc == 2);
    assert(cnet.group_sep == 1);
    assert((cnet.idx[0] == 3) && (cnet.idx[1] == 1) && (cnet.idx[2] == 2)
        && (cnet.idx[3] == 4) && (cnet.idx[4] == 5));
    assert((cnet.sep[0] == 0) && (cnet.sep[1] == 2) && (cnet.sep[2] == 5));
}

void test_squeeze2()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(3); 
    cnet.new_cable();
    cnet.append(1);
    cnet.new_cable();
    cnet.append(7);
    cnet.new_group();

    print_cablenet(cnet);
    int cnt = cnet.squeeze();
    print_cablenet(cnet);
    if (VERBOSE) printf("removed %d cables\n", cnt);
    assert(cnet.nc == 0);
    assert(cnet.group_sep == 0);
    assert(cnet.sep[0] == 0);
}

void test_squeeze3()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2);
    cnet.new_cable();
    cnet.append(1); cnet.append(2);
    cnet.new_group();
    cnet.new_cable();
    cnet.append(1); cnet.append(2);

    print_cablenet(cnet);
    assert(cnet.largest_index() == 2);
    int cnt = cnet.squeeze();
    print_cablenet(cnet);
    if (VERBOSE) printf("removed %d cables\n", cnt);
    assert(cnet.nc == 3);
    assert(cnet.group_sep == 2);
    assert(cnet.sep[3] == 6);
}

void test_remove_negative1()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2); cnet.append(-1); 
    cnet.new_cable();
    cnet.append(3); cnet.append(-1); cnet.append(4); cnet.append(5);
    cnet.new_group();
    cnet.new_cable();
    cnet.append(-1); cnet.append(6); cnet.append(7);
    print_cablenet(cnet);
    assert(cnet.largest_index() == 7);
    int cnt = cnet.remove_negative_indices();
    print_cablenet(cnet);
    if (VERBOSE) printf("removed %d nodes\n", cnt);
    assert(cnt == 3);
    assert(cnet.nc == 3);
    assert(cnet.group_sep == 2);
    assert((cnet.sep[0] == 0) && (cnet.sep[1] == 2) && (cnet.sep[2] == 5)
        && (cnet.sep[3] == 7));
    for (int i=0;i<7;i++) assert(cnet.idx[i] == i+1);
}

void test_remove_negative2()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(-1); cnet.append(-1); cnet.append(-1); 
    cnet.new_group();
    cnet.new_cable();
    cnet.append(1); cnet.append(-1); cnet.append(2);
    cnet.new_cable();
    cnet.append(-1);
    print_cablenet(cnet);
    assert(cnet.largest_index() == 2);
    int cnt = cnet.remove_negative_indices();
    print_cablenet(cnet);
    if (VERBOSE) printf("removed %d nodes\n", cnt);
    assert(cnt == 5);
    assert(cnet.nc == 3);
    assert(cnet.group_sep == 1);
    assert((cnet.sep[0] == 0) && (cnet.sep[1] == 0) && (cnet.sep[2] == 2)
        && (cnet.sep[3] == 2));
    for (int i=0;i<2;i++) assert(cnet.idx[i] == i+1);
}

void test_remove_negative3()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2); cnet.append(3); 
    cnet.new_group();
    print_cablenet(cnet);
    int cnt = cnet.remove_negative_indices();
    print_cablenet(cnet);
    if (VERBOSE) printf("removed %d nodes\n", cnt);
    assert(cnt == 0);
    assert(cnet.nc == 1);
    assert(cnet.group_sep == 1);
    assert((cnet.sep[0] == 0) && (cnet.sep[1] == 3) );
}

void test_remove_negative4()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(-1); 
    cnet.new_cable();
    cnet.append(-1); 
    print_cablenet(cnet);
    int cnt = cnet.remove_negative_indices();
    print_cablenet(cnet);
    if (VERBOSE) printf("removed %d nodes\n", cnt);
    assert(cnt == 2);
    assert(cnet.nc == 2);
    assert((cnet.sep[0] == 0) && (cnet.sep[1] == 0) );
}

void test_neighbors1()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2);  cnet.append(3); 
    cnet.new_cable();
    cnet.append(4); cnet.append(2);  cnet.append(5); 
    cnet.new_cable();
    cnet.append(5); cnet.append(6);
    print_cablenet(cnet);
    int nv = cnet.largest_index()+1;
    int* nn = new int [nv];
    cnet.count_neighbors(nv, nn);
    assert((nn[0] == 0) && (nn[1] == 1) && (nn[2] == 4) 
        && (nn[3] == 1) && (nn[4] == 1) && (nn[5] == 2) 
        && (nn[6] == 1));
    delete [] nn;
}

void test_neighbors2()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(0);
    cnet.new_cable();
    cnet.append(0); cnet.append(1);
    print_cablenet(cnet);
    int nv = cnet.largest_index()+1;
    int* nn = new int [nv];
    cnet.count_neighbors(nv, nn);
    assert((nn[0] == 2) && (nn[1] == 2));
    delete [] nn;
}

void test_neighbors3()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(0); cnet.append(2);
    cnet.new_cable();
    cnet.append(2); cnet.append(1);
    print_cablenet(cnet);
    int nv = cnet.largest_index()+1;
    int* nn = new int [nv];
    cnet.count_neighbors(nv, nn);
    assert((nn[0] == 2) && (nn[1] == 2) && (nn[2] == 2));
    delete [] nn;
}

void test_neighbors4()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(0); cnet.append(-1); cnet.append(1); 
    cnet.append(2); cnet.append(3);
    cnet.new_cable();
    cnet.append(2);
    print_cablenet(cnet);
    int nv = cnet.largest_index()+1;
    int* nn = new int [nv];
    cnet.count_neighbors(nv, nn);
    assert((nn[0] == 0) && (nn[1] == 1) && (nn[2] == 2) && (nn[3] == 1));
    delete [] nn;
}

void test_connectedness1()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(0); cnet.append(1); cnet.append(2); 
    cnet.new_cable();
    cnet.append(8);
    cnet.new_cable();
    cnet.append(2); cnet.append(3);
    cnet.new_cable();
    cnet.append(3);
    cnet.append(4);
    cnet.new_cable();
    cnet.append(6);
    cnet.append(5);
    print_cablenet(cnet);
    int ncomp;
    int cnt = cnet.remove_isolated_regions(ncomp);
    print_cablenet(cnet);
    assert(ncomp == 3);
    assert(cnet.size() == 7);
    assert(cnt == 3);
}

void test_connectedness2()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(0); cnet.append(1); cnet.append(2); 
    cnet.new_cable();
    cnet.append(4); cnet.append(5); cnet.append(6); cnet.append(7);
    print_cablenet(cnet);
    int ncomp;
    int cnt = cnet.remove_isolated_regions(ncomp);
    print_cablenet(cnet);
    assert(ncomp == 2);
    assert(cnet.size() == 4);
    assert(cnt == 3);
}

void test_connectedness3()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(4); cnet.append(5); cnet.append(6); cnet.append(7);
    cnet.new_cable();
    cnet.append(0); cnet.append(1); cnet.append(2); 
    print_cablenet(cnet);
    int ncomp;
    int cnt = cnet.remove_isolated_regions(ncomp);
    print_cablenet(cnet);
    assert(ncomp == 2);
    assert(cnet.size() == 4);
    assert(cnt == 3);
}

void test_connectedness4()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2); 
    cnet.new_cable();
    cnet.append(3); cnet.append(2);
    cnet.new_cable();
    cnet.append(4); cnet.append(1); cnet.append(8); 
    print_cablenet(cnet);
    int ncomp;
    int cnt = cnet.remove_isolated_regions(ncomp);
    print_cablenet(cnet);
    assert(ncomp == 1);
    assert(cnet.size() == 7);
    assert(cnt == 0);
}

void test_loose_ends1()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2); cnet.append(5); cnet.append(6); cnet.append(11); 
    cnet.new_cable();
    cnet.append(7); cnet.append(8); cnet.append(3); cnet.append(4);
    cnet.new_cable();
    cnet.append(1); cnet.append(3); cnet.append(9); 
    cnet.new_cable();
    cnet.append(10); cnet.append(2); cnet.append(4); 

    print_cablenet(cnet);
    int nv = cnet.largest_index()+1;
    int* nn = new int [nv];
    cnet.count_neighbors(nv, nn);
    int cnt = cnet.tag_loose_ends(nn);
    delete [] nn;
    print_cablenet(cnet);
    assert(cnt == 7);
}

void test_loose_ends2()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2); cnet.append(5); cnet.append(6); cnet.append(11); 
    cnet.new_cable();
    cnet.append(7); cnet.append(8); cnet.append(3); cnet.append(4);
    cnet.new_cable();
    cnet.append(1); cnet.append(3); cnet.append(9); 
    cnet.new_cable();
    cnet.append(10); cnet.append(2); cnet.append(4); 

    print_cablenet(cnet);
    int niter;
    int cnt = cnet.cut_loose_ends(niter);
    print_cablenet(cnet);
    assert(niter == 1);
    assert(cnt == 7);
}

void test_loose_ends3()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2); cnet.append(3); cnet.append(4);
    cnet.new_cable();
    cnet.append(2); cnet.append(4); cnet.append(5); cnet.append(6);
    cnet.new_cable();
    cnet.append(7); cnet.append(6); cnet.append(8);

    print_cablenet(cnet);
    int niter;
    int cnt = cnet.cut_loose_ends(niter);
    print_cablenet(cnet);
    assert(niter == 2);
    assert(cnet.size() == 5);
}

void test_loose_ends4()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2); cnet.append(3); cnet.append(4);
    cnet.new_cable();
    cnet.append(2); cnet.append(4); cnet.append(5); cnet.append(6);
    cnet.new_cable();
    cnet.append(7); cnet.append(5); cnet.append(8);

    print_cablenet(cnet);
    int niter;
    int cnt = cnet.cut_loose_ends(niter);
    print_cablenet(cnet);
    assert(niter == 2);
    assert(cnet.size() == 5);
}

void test_loose_ends5()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2); cnet.append(3); cnet.append(4);
    cnet.new_cable();
    cnet.append(2); cnet.append(4); cnet.append(5); cnet.append(6);
    cnet.new_cable();
    cnet.append(1); cnet.append(6); cnet.append(8);

    print_cablenet(cnet);
    int niter;
    int cnt = cnet.cut_loose_ends(niter);
    print_cablenet(cnet);
    assert(niter == 1);
    assert(cnet.size() == 10);
}

void test_loose_ends6()
{
    CableNetwork cnet;
    cnet.allocate(10, 20);
    cnet.clear();
    cnet.new_cable();
    cnet.append(1); cnet.append(2); cnet.append(3); cnet.append(4);
    cnet.new_cable();
    cnet.append(2); cnet.append(4); cnet.append(5); cnet.append(6);
    cnet.new_cable();
    cnet.append(7); cnet.append(6);
    cnet.new_cable();
    cnet.append(7); cnet.append(8);
    cnet.new_cable();
    cnet.append(9); cnet.append(8);

    print_cablenet(cnet);
    int niter;
    int cnt = cnet.cut_loose_ends(niter);
    print_cablenet(cnet);
    assert(niter == 4);
    assert(cnet.size() == 5);
    assert(cnt == 9);
}


int main()
{
    test_creation();
    test_squeeze0();
    test_squeeze0();
    test_squeeze1();
    test_squeeze2();
    test_squeeze3();
    test_remove_negative1();
    test_remove_negative2();
    test_remove_negative3();
    test_remove_negative4();
    test_neighbors1();
    test_neighbors2();
    test_neighbors3();
    test_neighbors4();
    test_connectedness1();
    test_connectedness2();
    test_connectedness3();
    test_connectedness4();
    test_loose_ends1();
    test_loose_ends2();
    test_loose_ends3();
    test_loose_ends4();
    test_loose_ends5();
    test_loose_ends6();
}