//
//  main.c
//  sfdp-lite
//
//  Created by Ciprian TEODOROV on 12/20/13.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>

#include "lib/mmio/mmio.h"

#include "graphviz/SparseMatrix.h"
#include "graphviz/spring_electrical.h"
#include "cairo/cairo.h"
#include "cairo/cairo-pdf.h"



#include "jet.h"

#define real double
#define FALSE 0
#define TRUE (!FALSE)

#define MAX(x,y) x<y?y:x
#define MIN(x,y) x<y?x:y

typedef struct pointf_s { double x, y; } pointf;

SparseMatrix makeSparseMatrix(int nz, int m, int n, int *i, int *j) {
    SparseMatrix theSparse;
    int theI;
    
    real *val = (real*)malloc((nz)*sizeof(real));
    if (val == NULL) {
        exit(3);
    }
    for (theI = 0; theI < nz; ++theI) {
        val[theI] = 1.0;
    }
    
    theSparse = SparseMatrix_from_coordinate_arrays(nz, m, n, i, j, val, MATRIX_TYPE_REAL, sizeof(real));
    
    if (val != NULL) free(val);
    
    return theSparse;
}


void draw_embedding(SparseMatrix A, real *x, const char *outFilePath){
    int i, j, *ia=A->ia, *ja = A->ja;
    real xsize, ysize, xmin, xmax, ymin, ymax, absxmin, absymin;
    real maximumDistance = 0;
    
    cairo_surface_t *surface;
    cairo_t *cr;
    
    int scale = 100;
    
    xmax = xmin = x[0];
    ymax = ymin = x[1];
    for (i = 0; i < A->m; i++){
        xmax = MAX(xmax, x[i*2]);
        xmin = MIN(xmin, x[i*2]);
        ymax = MAX(ymax, x[i*2+1]);
        ymin = MIN(ymin, x[i*2+1]);
    }
    xsize = xmax-xmin;
    ysize = ymax-ymin;
    xsize = MAX(xsize, ysize);
    absxmin = (xmin < 0 ? -xmin : xmin)*scale;
    absymin = (ymin < 0 ? -ymin : ymin)*scale;
    
    //find the maximum distance
    for (i = 0; i < A->m; i++){
        for (j = ia[i]; j < ia[i+1]; j++){
            if (ja[j] == i) continue;
            real dist = fabs(x[i*2+0] - x[ja[j]*2+0]) + fabs(x[i*2+1] - x[ja[j]*2+1]);
            
            maximumDistance = MAX(dist, maximumDistance);
        }
    }

    printf("drawing size (%lf, %lf)", xsize*scale, ysize*scale);
    
    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, xsize*scale, ysize*scale);
    //cairo_pdf_surface_create(outFilePath, xsize*scale, ysize*scale);
    cr = cairo_create (surface);
    
    
    
    for (i = 0; i < A->m; i++){
        for (j = ia[i]; j < ia[i+1]; j++){
            if (ja[j] == i) continue;
            
            real dist = fabs(x[i*2+0] - x[ja[j]*2+0]) + fabs(x[i*2+1] - x[ja[j]*2+1]);
            //normalize dist to [0,1] interval
            dist = dist / maximumDistance;
            //get an integer id mapping for dist
            int id = floor(dist*1024);
            //get the correct id from jet which is in reverse order
            id = 1024 - id;
            
            cairo_set_source_rgba (cr, jet[id][0], jet[id][1], jet[id][2], 0.5);

            cairo_move_to(cr, x[i*2+0]*scale + absxmin, x[i*2+1]*scale + absymin);
            cairo_line_to(cr, x[ja[j]*2+0]*scale + absxmin, x[ja[j]*2+1]*scale + absymin);
            
            cairo_stroke(cr);
        }
    }

//    cairo_set_source_rgb(cr, 0, 0, 1);
//    cairo_arc(cr, x[0]*scale*absxmin, x[1]*scale + absymin, 100, 0, 2*M_PI);
//    cairo_fill(cr);

    //cairo_surface_finish(surface);
    cairo_surface_write_to_png(surface , outFilePath);
}

#define PDF_OBP "./examples/obp.pdf"
#define PDF_OBP1 "./examples/obp1.pdf"
#define PDF_OBP2 "./examples/obp2.png"
#define PDF_OBP3 "./examples/obp3.png"
#define PDF_OBP4 "./examples/obp4.png"
#define PDF_UTM "./examples/utm1700b.pdf"
#define PDF_FILE_PATH PDF_OBP4
#define DOLAYOUT FALSE
#define SPANNING_TREE FALSE

void doSfdpLayout(SparseMatrix matrix, spring_electrical_control ctrl, short doLayout, real *positions, FILE *outPositions) {
    int flag;
    clock_t runtime;
    
    if (doLayout) {
        //start layouting
        printf("starting SFDP layout\n");
        runtime = clock();
        multilevel_spring_electrical_embedding(2, matrix, NULL, ctrl, NULL, NULL, positions, 0, NULL, &flag);
        runtime = clock() - runtime;
        printf("Finished layout in %f seconds\n", (float)runtime / CLOCKS_PER_SEC);
        
        //write position in the file
        runtime = clock();
        int i;
        for (i = 0; i < matrix->n; ++i) {
            real *npos = positions + 2 * i;
            fprintf(outPositions, "%lg %lg\n", npos[0], npos[1]);
        }
        runtime = clock() - runtime;
        printf("Wrote positions took %f seconds\n", (float)runtime / CLOCKS_PER_SEC);
    }
    
    //draw the embedding with cairo
    printf("starting drawing\n");
    runtime = clock();
    draw_embedding(matrix, positions, PDF_FILE_PATH);
    runtime = clock() - runtime;
    printf("Drawing layout took %f seconds\n", (float)runtime / CLOCKS_PER_SEC);
    
    //free the matrix
    SparseMatrix_delete(matrix);
}

static void tuneControl(spring_electrical_control ctrl) {
    ctrl->random_seed = (unsigned) getpid() ^ (unsigned) time(NULL);
    ctrl->K = -1;
    ctrl->p = -1.0 * -AUTOP;
    ctrl->multilevels = INT_MAX;
    ctrl->smoothing = SMOOTHING_NONE;
    ctrl->tscheme = QUAD_TREE_NORMAL;
    ctrl->method = METHOD_SPRING_ELECTRICAL;
    ctrl->beautify_leaves = FALSE;
    ctrl->do_shrinking = TRUE;
    ctrl->rotation = 0.0;
    ctrl->edge_labeling_scheme = 0;
    ctrl->overlap = -1; //turn off overlap removal
}

typedef struct cell_s cell_t;
struct cell_s {
    int node;
    cell_t *next;
};

typedef struct list_s list_t;
struct list_s {
    cell_t *first;
    cell_t *last;
};

int dequeue(list_t *inList) {
    if (inList == NULL || inList->first == NULL) {
        return -1;
    }
    cell_t *cell = inList->first;
    
    if (inList->first == inList->last) {
        inList->first = inList->last = NULL;
    }
    else {
        inList->first = cell->next;
    }
    int value = cell->node;
    
    free(cell);
    return value;
}
#include <assert.h>
list_t* enqueue(list_t *inList, int value) {
    
    if (inList == NULL) {
        list_t *theList = malloc(sizeof(list_t));
        cell_t *cell = (cell_t*)malloc(sizeof(cell_t));
        cell->node = value;
        cell->next = NULL;
        theList->first = theList->last = cell;
        return theList;
    }
    
    cell_t *cell = (cell_t*)malloc(sizeof(cell_t));
    
    cell->node = value;
    cell->next = NULL;
    
    if (inList->last == NULL) {
        assert(inList->first == NULL);
        inList->first = inList->last = cell;
        return inList;
    }
    
    inList->last->next = cell;
    inList->last = cell;
    return inList;
}

SparseMatrix computeSpanningTree(SparseMatrix matrix) {
    SparseMatrix spanning;
    int j, *ia=matrix->ia, *ja = matrix->ja;
    
    char *known;
    list_t *theQueue;
    
    int *adj[2];
    adj[0] = (int*)malloc(matrix->nz*sizeof(int));
    adj[1] = (int*)malloc(matrix->nz*sizeof(int));

    int nedge = 0;
    
    known = (char *)calloc(matrix->n, sizeof(char));
    
    theQueue = enqueue(NULL, 0);
    
    while (theQueue->first != NULL) {
        int current = dequeue(theQueue);
        
        for (j = ia[current]; j < ia[current+1]; j++){
            int to = ja[j];
            //ignore self links
            if (to == current) continue;
            
            if (!known[to]) {
                known[to] = 1;
                theQueue = enqueue(theQueue, to);
                
                adj[0][nedge] = current;
                adj[1][nedge++] = to;
            }
        }
    }
    
    free(theQueue);
    free(known);
    
    spanning = makeSparseMatrix(nedge, matrix->m, matrix->n, adj[0], adj[1]);
    free(adj[0]);
    free(adj[1]);
    
    
    return spanning;
}

int main(int argc, const char * argv[])
{

    MM_typecode theMatcode;
    int         theM,       theN, theNZ;
    int        *theI,      *theJ;
    double     *theValues;
    FILE *thePositionsFile = NULL;
    SparseMatrix theSparse;
    spring_electrical_control ctrl = spring_electrical_control_new();
    
    if (argc < 2 || (argc > 2 && argc < 4)) {
        fprintf(stderr, "Usage %s [matrix-market-filename] -o [positions-file]\n", argv[0]);
        exit(1);
    }
    //read the graph from matrix market file
    if (mm_read_mtx_crd(argv[1],
                    &theM, &theN, &theNZ,
                    &theI, &theJ, &theValues,
                        &theMatcode) != 0) {
        fprintf(stderr, "Failed reading file %s\n", argv[1]);
        exit(1);
    }
    
    printf("Finished reading %d by %d matrix and %d nonzeros\n", theM, theN, theNZ);
    
    //make sparse matrix
    theSparse = makeSparseMatrix(theNZ, theM, theN, theI, theJ);
    
    if (SPANNING_TREE) {
        SparseMatrix spanning;
        spanning = computeSpanningTree(theSparse);
        SparseMatrix_delete(theSparse);
        theSparse = spanning;
    }
    
    if (theI != NULL) free(theI);
    if (theJ != NULL) free(theJ);
    if (theValues != NULL) free(theValues);
    
    printf("Sparse created\n");
    


    
    tuneControl(ctrl);
    //only one connected component in graphs from OBP
    
    //initial node position is (0,0)
    real *positions = (real*)calloc(2 * theN, sizeof(real));
    if (DOLAYOUT) {
        //prepare the output file
        thePositionsFile = argc == 2 ? stdout : fopen(argv[3], "w");
        
        if (thePositionsFile != stdout && thePositionsFile == NULL) {
            fprintf(stderr, "Could not open the specified output file\n");
            exit(1);
        }

    }
    else {
        FILE *f = fopen(argv[3], "r");
        
        for (int i = 0; i < theSparse->n; ++i) {
            real *npos = positions + 2 * i;
            fscanf(f, "%lg %lg", &npos[0], &npos[1]);
            //printf("%lg %lg\n", npos[0], npos[1]);
        }
        
        fclose(f);
    }
    doSfdpLayout(theSparse, ctrl, DOLAYOUT, positions, thePositionsFile);
    
    spring_electrical_control_delete(ctrl);
    
    if (thePositionsFile != stdout) fclose(thePositionsFile);
    
    free(positions);
    printf("Finished saving\n");

    return 0;
}


