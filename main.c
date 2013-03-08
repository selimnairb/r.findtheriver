
/****************************************************************************
 *
 * MODULE:       r.findtheriver
 * AUTHOR(S):    Brian Miles - brian_miles@unc.edu
 *               with hints from: Glynn Clements - glynn gclements.plus.com
 * PURPOSE:      Finds the nearest stream pixel to a coordinate pair using
 * 				 an upstream accumulating area map.
 *
 * COPYRIGHT:    (C) 2013 by the University of North Carolina at Chapel Hill
 *
 *               This program is free software under the GNU General Public
 *   	    	 License (>=v2). Read the file COPYING that comes with GRASS
 *   	    	 for details.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <grass/gis.h>
#include <grass/glocale.h>

#include "point_list.h"

#define THRESHOLD_DISTANCE 100
#define DEFAULT_COORD_SEP ' '

/* 
 * global function declaration 
 */
PointList_t *findStreamPixelsInWindow(int fd, RASTER_MAP_TYPE dataType, int windowSize, int threshold,
		int nrows_less_one, int ncols_less_one,
		int currRow, int currCol);

/*
 * function definitions 
 */
PointList_t *findStreamPixelsInWindow(int fd, RASTER_MAP_TYPE dataType, int windowSize, int threshold,
		int nrows_less_one, int ncols_less_one,
		int currRow, int currCol) {
	assert( DCELL_TYPE == dataType );
	assert( windowSize & 1 );
	assert( threshold > 0 );

	PointList_t *streamPixels = NULL;
	DCELL centralValue, tmpValue;
	double logCentralValue, logTmpValue;
	void *tmpRow = G_allocate_raster_buf(dataType);

	// Get value of central cell
	if (G_get_raster_row(fd, tmpRow, currRow, dataType) < 0) {
		G_fatal_error(_("Unable to read raster row %d"), currRow);
	}

	switch (dataType) {
	case CELL_TYPE:
		centralValue = (double)((CELL *) tmpRow)[currCol];
		break;
	case FCELL_TYPE:
		centralValue = (double)((FCELL *) tmpRow)[currCol];
		break;
	case DCELL_TYPE:
		centralValue = (double)((DCELL *) tmpRow)[currCol];
		break;
	}
	logCentralValue = log10(centralValue);

	// Define window bounds
	int windowOffset = ( windowSize - 1 ) / 2;
	int minCol = currCol - windowOffset;
	if ( minCol < 0 ) minCol = 0;
	int maxCol = currCol + windowOffset;
	if ( maxCol > ncols_less_one ) maxCol = ncols_less_one;
	int minRow = currRow - windowOffset;
	if ( minRow < 0 ) minRow = 0;
	int maxRow = currRow + windowOffset;
	if ( maxRow > nrows_less_one ) maxRow = nrows_less_one;

	// Search for stream pixels within the window
	int row, col;
	for ( row = minRow ; row <= maxRow ; row++ ) {
		// Get the current row
		if (G_get_raster_row(fd, tmpRow, row, dataType) < 0) {
				G_fatal_error(_("Unable to read raster row %d"), row);
		}
		for ( col = minCol ; col <= maxCol ; col++ ) {
			switch (dataType) {
			case CELL_TYPE:
				tmpValue = (double)((CELL *) tmpRow)[currCol];
				break;
			case FCELL_TYPE:
				tmpValue = (double)((FCELL *) tmpRow)[currCol];
				break;
			case DCELL_TYPE:
				tmpValue = (double)((DCELL *) tmpRow)[currCol];
				break;
			}
			logTmpValue = log10(tmpValue);
			// Test for nearby pixels that are stream pixels when compared to the central pixel
			//fprintf(stderr, "logTmpValue: %f, logCentralValue: %f\n",
			//		logTmpValue, logCentralValue);
			if ( (logTmpValue - logCentralValue) > threshold ) {
				// Add to list of stream pixels
				if ( NULL == streamPixels ) {
					streamPixels = createList(col, row);
				} else {
					appendPoint(streamPixels, col, row);
				}
			}
		}
	}
	G_free(tmpRow);
	return streamPixels;
}

int main(int argc, char *argv[])
{
	struct Cell_head cellhd;	/* it stores region information,
				   and header information of rasters */
	char name[GNAME_MAX];			/* input raster name */
	char *mapset;		/* mapset name */
	int nrows, ncols;
	int rowIdx, colIdx, nrows_less_one, ncols_less_one, total;
	int infd;		/* file descriptor */
	int quiet;
	RASTER_MAP_TYPE data_type;	/* type of the map (CELL/DCELL/...) */
	struct GModule *module;	/* GRASS module for parsing arguments */

	struct Flag *flag1;		/* flags */
	struct Option *optInput, *optWindow, *optThreshold, *optE, *optN, *optSep;
	double E, N;
	struct Cell_head window;
	char *buff;
	int windowSize, threshold;
	size_t SEP_SIZE = 1;
	char sep[SEP_SIZE];

	/* initialize GIS environment */
	G_gisinit(argv[0]);		/* reads grass env, stores program name to G_program_name() */

	/* initialize module */
	module = G_define_module();
	module->keywords = _("raster, keyword2, keyword3");
	module->description = _("Find the stream pixel nearest the input coordinate");

	/* Define command options */
	optInput = G_define_option();
	optInput->key = "accumulation";
	optInput->type = TYPE_STRING;
	optInput->required = YES;
	optInput->gisprompt = "old,cell,raster";
	optInput->description = _("Name of input upstream accumulation area raster map");

	optWindow = G_define_option();
	optWindow->key = "window";
	optWindow->type = TYPE_INTEGER;
	optWindow->key_desc = "x";
	optWindow->multiple = NO;
	optWindow->required = NO;
	optWindow->description = _("The size of the window, in pixels, to search in for stream pixels.  Must be an odd integer.  If not supplied, window will be inferred based on raster resolution.");

	optThreshold = G_define_option();
	optThreshold->key = "threshold";
	optThreshold->type = TYPE_INTEGER;
	optThreshold->key_desc = "x";
	optThreshold->multiple = NO;
	optThreshold->required = NO;
	optThreshold->description = _("The threshold for distinguishing log(UAA) values of stream and non-stream pixels.  If not supplied, threshold will be inferred from minimum and maximum raster values.");

	optE = G_define_option();
	optE->key = "easting";
	optE->type = TYPE_STRING;
	optE->key_desc = "x";
	optE->multiple = NO;
	optE->required = YES;
	optE->description = _("The map E grid coordinates");

	optN = G_define_option();
	optN->key = "northing";
	optN->type = TYPE_STRING;
	optN->key_desc = "y";
	optN->multiple = NO;
	optN->required = YES;
	optN->description = _("The map N grid coordinates");

	optSep = G_define_option();
	optSep->key = "separator";
	optSep->type = TYPE_STRING;
	optSep->key_desc = "y";
	optSep->multiple = NO;
	optSep->required = NO;
	optSep->description = _("Coordinate separator. Defaults to ' '. Must be 1 character in length.");

	/* Define the different flags */
	flag1 = G_define_flag();
	flag1->key = 'q';
	flag1->description = _("Quiet");

	/* options and flags parser */
	if (G_parser(argc, argv))
		exit(EXIT_FAILURE);

	if (G_get_window(&window) < 0) {
		G_asprintf(&buff, _("Unable to read current window parameters"));
		G_fatal_error(buff);
	}

	/* stores options and flags to variables */
	strncpy(name, optInput->answer, GNAME_MAX);
	quiet = (flag1->answer);

	/* returns NULL if the map was not found in any mapset,
	 * mapset name otherwise */
	mapset = G_find_cell2(name, "");
	if (mapset == NULL) {
		G_fatal_error(_("Raster map <%s> not found"), name);
	}

	/* Get raster metadata */
	if (G_get_cellhd(name, mapset, &cellhd) < 0) {
		G_fatal_error(_("Unable to read file header of <%s>"), name);
	}

	if ( NULL != optWindow->answer ) {
		windowSize = atoi(optWindow->answer);
	} else {
		// Determine window size
		double cellRes = (cellhd.ew_res + cellhd.ns_res) / 2;
		windowSize = THRESHOLD_DISTANCE / cellRes;
		if ( !(windowSize & 1)  ) windowSize++;
	}
	if ( (windowSize < 2) || !(windowSize & 1) ) {
		G_warning(_("Invalid window size %s.  Window size must be an odd integer >= 3\n"), optWindow->answer);
		G_usage();
		exit(EXIT_FAILURE);
	}
	if ( !quiet ) {
		fprintf(stderr, "Stream search window size %d\n", windowSize);
	}
	if ( NULL != optThreshold->answer ) {
		threshold = atoi(optThreshold->answer);
	} else {
		// Determine threshold
		struct FPRange *range = (struct FPRange *)malloc(sizeof(struct FPRange));
		if ( G_read_fp_range(name, mapset, range) < 0 ) {
			G_fatal_error(_("Unable to determine range of raster map <%s>"), name);
		}
		double min = range->min;
		double max = range->max;
		free(range); // eggs
		if ( min <= 0 ) min = 1;
		double logMin = log10(min);
		double logMax = log10(max);
		threshold = (logMax - logMin) - 1;
	}
	if ( threshold < 1 ) {
		G_warning(_("Invalid threshold %s.  Window size must be an integer >= 1\n"), optThreshold->answer);
		G_usage();
		exit(EXIT_FAILURE);
	}
	if ( !quiet ) {
		fprintf(stderr, "Stream log-difference threshold %d\n", threshold);
	}

	if (!G_scan_easting(*optE->answers, &E, G_projection())) {
		G_warning(_("Illegal east coordinate <%s>\n"), optE->answer);
		G_usage();
		exit(EXIT_FAILURE);
	}
	if (!G_scan_northing(*optN->answers, &N, G_projection())) {
		G_warning(_("Illegal north coordinate <%s>\n"), optN->answer);
		G_usage();
		exit(EXIT_FAILURE);
	}
	if ( !quiet ) {
		fprintf(stderr, "Input coordinates, easting %f, northing %f\n", E, N);
	}

	if ( NULL == optSep->answer) {
		sep[0] = DEFAULT_COORD_SEP;
	} else {
		strncpy(&sep, optSep->answer, SEP_SIZE);
	}

	/* Determine the inputmap type (CELL/FCELL/DCELL) */
	data_type = G_raster_map_type(name, mapset);

	/* Open the raster - returns file destriptor (>0) */
	if ((infd = G_open_cell_old(name, mapset)) < 0)
		G_fatal_error(_("Unable to open raster map <%s>"), name);

	G_get_set_window(&window);
	nrows = G_window_rows();
	ncols = G_window_cols();
	total = nrows * ncols;
	nrows_less_one = nrows - 1;
	ncols_less_one = ncols - 1;

	rowIdx = (int)G_northing_to_row(N, &window);
	colIdx = (int)G_easting_to_col(E, &window);

	double currNearestE, prevNearestE, currNearestN, prevNearestN;
	PointList_t *streamPixels = findStreamPixelsInWindow(infd, data_type, windowSize, threshold,
			nrows_less_one, ncols_less_one,
			rowIdx, colIdx);
	if ( !quiet ) {
		fprintf(stderr, "Stream pixels: ");
		printList(stderr, streamPixels, " ");
		fprintf(stderr, "\n");
	}
	PointList_t *nearestStreamPixel = findNearestPoint(streamPixels, colIdx, rowIdx);

	if ( NULL != nearestStreamPixel ) {
		double nearestEasting = G_col_to_easting(nearestStreamPixel->col, &window);
		double nearestNorthing = G_row_to_northing(nearestStreamPixel->row, &window);
		printf("%f%s%f\n", nearestEasting, sep, nearestNorthing);
	}

	// Clean up
	destroyList(streamPixels);

	/* closing raster maps */
	G_close_cell(infd);

	exit(EXIT_SUCCESS);
}
