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

#ifndef boolean
#define boolean _Bool
#endif

#define real double

typedef struct {
    char   *inputFile;
    char   *outputFile;
    char   *positionFile;
    char   *traceFile;
    char   *layoutFile;
    char   *initialPositionsFile;
    int     scaleX;
    int     scaleY;
    float   alpha;
    boolean isJet;
    boolean isSpanning;
} arguments_t;

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

#define JET_COLORS
//#define SHOW_TRACE

#ifdef SHOW_TRACE
//#define TRACE_SIZE 40
//int trace[40][2] = {{258, 525},{80, 345},{895, 347},{599, 616},{624, 643},{480, 389},{921, 555},{555, 561},{822, 67},{58, 318},{997, 242},{879, 894},{135, 152},{981, 997},{67, 80},{900, 921},{643, 900},{189, 480},{242, 258},{894, 135},{700, 718},{238, 879},{561, 636},{525, 540},{636, 189},{540, 810},{41, 58},{152, 416},{318, 335},{335, 599},{416, 432},{347, 41},{389, 684},{345, 359},{684, 238},{810, 822},{616, 879},{359, 624},{432, 700},{718, 981}};

#define TRACE_SIZE 66
int trace[66][2] = {{1, 2},{2, 5},{5, 10},{10, 13},{13, 16},{16, 19},{19, 22},{22, 25},{25, 28},{28, 31},{31, 34},{34, 37},{37, 40},{40, 43},{43, 46},{46, 49},{49, 52},{52, 55},{55, 58},{58, 61},{61, 64},{64, 67},{67, 70},{70, 73},{73, 76},{76, 79},{79, 82},{82, 85},{85, 88},{88, 91},{91, 94},{94, 97},{97, 102},{102, 112},{112, 136},{136, 160},{160, 190},{190, 226},{226, 268},{268, 316},{316, 370},{370, 397},{397, 424},{424, 451},{451, 478},{478, 505},{505, 532},{532, 559},{559, 586},{586, 613},{613, 640},{640, 667},{667, 694},{694, 721},{721, 748},{748, 775},{775, 802},{802, 829},{829, 856},{856, 883},{883, 910},{910, 937},{937, 964},{964, 991},{991, 1022},{1022, 112}};
#endif

int * read_trace(const char *inTraceFile, int *outTraceSize) {
    int * theTrace;
    FILE *f = fopen(inTraceFile, "r");
    fscanf(f, "%d", outTraceSize);
    theTrace = (int*)calloc(*outTraceSize*2, sizeof(real));
    for (int i = 0; i < *outTraceSize; ++i) {
        int *npos = theTrace + 2 * i;
        fscanf(f, "%d %d", &npos[0], &npos[1]);
    }
    
    fclose(f);
    return theTrace;
}

void draw_embedding(SparseMatrix A, real *x, arguments_t *args){
    int i, j, *ia=A->ia, *ja = A->ja;
    real xsize, ysize, xmin, xmax, ymin, ymax, absxmin, absymin;
    real maximumDistance = 0;
    
    cairo_surface_t *surface;
    cairo_t *cr;
    
    int scalex = args->scaleX;
    int scaley = args->scaleY;
    
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
    absxmin = (xmin < 0 ? -xmin : xmin)*scalex;
    absymin = (ymin < 0 ? -ymin : ymin)*scaley;
    
    //find the maximum distance
    for (i = 0; i < A->m; i++){
        for (j = ia[i]; j < ia[i+1]; j++){
            if (ja[j] == i) continue;
            real dist = fabs(x[i*2+0] - x[ja[j]*2+0]) + fabs(x[i*2+1] - x[ja[j]*2+1]);
            
            maximumDistance = MAX(dist, maximumDistance);
        }
    }

    printf("drawing size (%lf, %lf)", xsize*scalex, ysize*scaley);
    
    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, xsize*scalex, ysize*scaley);
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
     
            if (args->isJet) {
                cairo_set_source_rgba (cr, jet[id][0], jet[id][1], jet[id][2], args->alpha);
            } else {
                cairo_set_source_rgba (cr, 0.5, 0.5, 0.5, args->alpha);
            }
            cairo_set_line_width(cr,0.8);
            cairo_move_to(cr, x[i*2+0]*scalex + absxmin, x[i*2+1]*scaley + absymin);
            cairo_line_to(cr, x[ja[j]*2+0]*scalex + absxmin, x[ja[j]*2+1]*scaley + absymin);
            
            cairo_stroke(cr);
        }
    }

    if (args->traceFile == NULL) {
        cairo_set_source_rgb(cr, 0, 1, 0);
        cairo_arc(cr, x[0]*scalex+absxmin, x[1]*scaley + absymin, 0.05*scaley, 0, 2*M_PI);
        cairo_fill(cr);
    } else {
        int   theTraceSize;
        int  *theTrace      = read_trace(args->traceFile, &theTraceSize);
        for (i = 0; i < theTraceSize; i++) {
            int *npos = theTrace + 2 * i;
            cairo_set_source_rgb(cr, 0, 0, 0);
            cairo_set_line_width(cr,3.0);
            cairo_move_to(cr, x[(npos[0]-1)*2+0]*scalex + absxmin, x[(npos[0]-1)*2+1]*scaley + absymin);
            cairo_line_to(cr, x[(npos[1]-1)*2+0]*scalex + absxmin, x[(npos[1]-1)*2+1]*scaley + absymin);
        
            cairo_stroke(cr);
        }
    }
    //cairo_surface_finish(surface);
    cairo_surface_write_to_png(surface , args->outputFile);
}

void readPositions(const char *inFile, real *positions, int size) {
    //the layout is precomputed read it from the file
    FILE *f = fopen(inFile, "r");
    
    for (int i = 0; i < size; ++i) {
        real *npos = positions + 2 * i;
        fscanf(f, "%lg %lg", &npos[0], &npos[1]);
    }
    
    fclose(f);
}

void writePositions(const char *outFile, real *positions, int size) {
    clock_t runtime;
    //prepare the output file
    FILE *theFile = fopen(outFile, "w");
    
    if (theFile == NULL) {
        fprintf(stderr, "Could not write positions in %s file\n", outFile);
        exit(1);
    }

    //write position in the file
    runtime = clock();
    int i;
    for (i = 0; i < size; ++i) {
        real *npos = positions + 2 * i;
        fprintf(theFile, "%lg %lg\n", npos[0], npos[1]);
    }
    runtime = clock() - runtime;
    printf("Wrote positions took %f seconds\n", (float)runtime / CLOCKS_PER_SEC);
    if (theFile != NULL) fclose(theFile);
}

void doSfdpLayout(SparseMatrix matrix,
                  spring_electrical_control ctrl,
                  arguments_t *args) {
    int flag;
    clock_t runtime;
    
    //initial node position is (0,0)
    real *positions = (real*)calloc(2 * matrix->n, sizeof(real));
    if (args->layoutFile == NULL) {
        if (args->initialPositionsFile != NULL) {
            readPositions(args->initialPositionsFile, positions, matrix->n);
        }
        
        //start layouting
        printf("starting SFDP layout\n");
        runtime = clock();
        multilevel_spring_electrical_embedding(2, matrix, NULL, ctrl, NULL, NULL, positions, 0, NULL, &flag);
        runtime = clock() - runtime;
        printf("Finished layout in %f seconds\n", (float)runtime / CLOCKS_PER_SEC);
       
        if (args->positionFile) {
            writePositions(args->positionFile, positions, matrix->n);
        } 
       
    } else {
        //the layout is precomputed read it from the file
        readPositions(args->layoutFile, positions, matrix->n);
    }
    
    //draw the embedding with cairo
    printf("starting drawing\n");
    runtime = clock();
    draw_embedding(matrix, positions, args);
    runtime = clock() - runtime;
    printf("Drawing layout took %f seconds\n", (float)runtime / CLOCKS_PER_SEC);
    
    //free the matrix
    SparseMatrix_delete(matrix);
    
    free(positions);
}

static void tuneControl(spring_electrical_control ctrl) {
    ctrl->random_seed = (unsigned) getpid() ^ (unsigned) time(NULL);
    ctrl->K = -1;
    ctrl->q = -2.0;
    ctrl->p = -1.0 * -AUTOP;
    //ctrl->tol = 0.00000001;
    ctrl->multilevels = INT_MAX;
    ctrl->smoothing = SMOOTHING_NONE;
    ctrl->tscheme = QUAD_TREE_NORMAL;
    ctrl->method = METHOD_SPRING_ELECTRICAL;
    ctrl->beautify_leaves = TRUE;
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

int sfdp_main(arguments_t *args)
{
    MM_typecode theMatcode;
    int         theM,       theN, theNZ;
    int        *theI,      *theJ;
    double     *theValues;
    SparseMatrix theSparse;
    spring_electrical_control ctrl = spring_electrical_control_new();
    
    //read the graph from matrix market file
    if (mm_read_mtx_crd(args->inputFile,
                        &theM, &theN, &theNZ,
                        &theI, &theJ, &theValues,
                        &theMatcode) != 0) {
        fprintf(stderr, "Failed reading file %s\n", args->inputFile);
        exit(1);
    }
    
    printf("Finished reading %d by %d matrix and %d nonzeros\n", theM, theN, theNZ);
    
    //make sparse matrix
    theSparse = makeSparseMatrix(theNZ, theM, theN, theI, theJ);
    
    if (args->isSpanning) {
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
    
    doSfdpLayout(theSparse, ctrl, args);
    
    spring_electrical_control_delete(ctrl);
    
    return 0;
}



void
usage(int estatus) {
    printf("Usage: sfdp-lite [-options] <mtx-file>\n");
    printf(" where the options include:\n");
    printf("\t[-t file | --show-trace file]\t specify a trace file\n");
    printf("\t[-o file | --output file]\t specify the output filename [default=output.png]\n");
    printf("\t[-p file | --position file]\t specify the position filename\n");
    printf("\t[-r file | --read-layout file]\t read a layout from the file");
    printf("\t[-i file | --initial-positions file]\t read initial positions from file");
    printf("\t[-x int | --scale-x int]\t\t scale x\n");
    printf("\t[-y int | --scale-y int]\t\t scale y\n");
    printf("\t[-a float | --alpha float]\t\t alpha in [0,1]\n");
    printf("\t[-j | --jet]\t\t\t\t\t\t jet colors\n");
    printf("\t[-s | --spanning]\t\t\t\t compute spanning tree layout\n");
    printf("\t[-v | --version]\t\t\t\t\t\t print sfpd-lite version and exit\n");
    printf("\t[-h | --help]\t\t\t\t\t\t print this help\n");
    
    exit(estatus);
}
#define VERSION_NUM "v1.0"
void
print_version() {
    printf("sfdp-lite %s\n", VERSION_NUM);
}

int main(int argc, const char * argv[]) {
    int i;
    arguments_t args = {NULL, NULL, NULL, NULL, NULL, NULL, 100, 100, 0.5, TRUE, FALSE};
    for (i=1; i<argc-1;i++){
        if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i],"--output") == 0) {
            if (i+1 > argc) { usage(1); }
            args.outputFile = (char*)argv[i+1];
            i++;
        } else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i],"--show-trace") == 0) {
            if (i+1 > argc) { usage(1); }
            args.traceFile = (char*)argv[i+1];
            i++;
        } else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i],"--read-layout") == 0) {
            if (i+1 > argc) { usage(1); }
            args.layoutFile = (char*)argv[i+1];
            i++;
        } else if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i],"--initial-positions") == 0) {
            if (i+1 > argc) { usage(1); }
            args.initialPositionsFile = (char*)argv[i+1];
            i++;
        } else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i],"--position") == 0) {
            if (i+1 > argc) { usage(1); }
            args.positionFile = (char*)argv[i+1];
            i++;
        } else if (strcmp(argv[i], "-x") == 0 || strcmp(argv[i],"--scale-x") == 0) {
            if (i+1 > argc) { usage(1); }
            args.scaleX = atoi(argv[i+1]);
            i++;
        } else if (strcmp(argv[i], "-y") == 0 || strcmp(argv[i],"--scale-y") == 0) {
            if (i+1 > argc) { usage(1); }
            args.scaleY = atoi(argv[i+1]);
            i++;
        } else if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i],"--alpha") == 0) {
            if (i+1 > argc) { usage(1); }
            args.alpha = atoi(argv[i+1]);
            i++;
        } else if (strcmp(argv[i], "-j") == 0 || strcmp(argv[i],"--no-jet") == 0) {
            if (i+1 > argc) { usage(1); }
            args.isJet = FALSE;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i],"--spanning") == 0) {
            if (i+1 > argc) { usage(1); }
            args.isSpanning = TRUE;
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i],"--version") == 0) {
            print_version();
            exit(0);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i],"--help") == 0) {
            usage(0);
        }
    }
    args.inputFile = (char*) argv[i];
    
    return sfdp_main(&args);
}



