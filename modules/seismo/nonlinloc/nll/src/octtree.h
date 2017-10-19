/*
 * Copyright (C) 1999-2005 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


/*  octree.h

	include file for octree search

*/



/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#ifdef EXTERN_MODE
#define	EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif

	/* misc defines */

#ifndef SMALL_DOUBLE
#define SMALL_DOUBLE 1.0e-20
#endif
#ifndef LARGE_DOUBLE
#define LARGE_DOUBLE 1.0e20
#endif
#ifndef VERY_SMALL_DOUBLE
#define VERY_SMALL_DOUBLE 1.0e-30
#endif
#ifndef VERY_LARGE_DOUBLE
#define VERY_LARGE_DOUBLE 1.0e30
#endif


/*------------------------------------------------------------/ */
/* structures */
/*------------------------------------------------------------/ */


/* octree node */

typedef struct octnode* OctNodePtr;
typedef struct octnode
{
	OctNodePtr parent;		/* parent node */
	Vect3D center;			/* absolute coordinates of center */
	Vect3D ds;			/* length of sides */
	double value;			/* node value */
	OctNodePtr child[2][2][2];	/* child nodes */
	char isLeaf;			/* leaf flag, 1=leaf */
	void *pdata;		/* additional data */
} OctNode;



/* 3D tree with Nx, Ny, Nz arbitrary */

typedef struct
{
	OctNode**** nodeArray;		/* parent nodes */
	int data_code;		/* data type code, application dependent */
	int numx, numy, numz;		/* grid size */
 	Vect3D orig; 	/* orig (km) */
	Vect3D ds;		/* len side (km) */
	double integral;
}
Tree3D;


/* structure for storing results */

typedef struct resultTreeNode* ResultTreeNodePtr;
typedef struct resultTreeNode {
	ResultTreeNodePtr left;		/* address of left node */
	ResultTreeNodePtr right;	/* address of right node */
	double value;			/* sort value */
	double volume;		/* volume, node volume depends on geometry in physical space, may not be dx*dy*dz */
	OctNode* pnode;			/* correspnding octree node */
} ResultTreeNode;



/* */
/*------------------------------------------------------------/ */



/*------------------------------------------------------------/ */
/* globals  */
/*------------------------------------------------------------/ */

//EXTERN_TXT char fn_control[MAXLINE];	/* control file name */

/* */
/*------------------------------------------------------------/ */




/*------------------------------------------------------------/ */
/* function declarations */
/*------------------------------------------------------------/ */

Tree3D* newTree3D(int data_code, int numx, int numy, int numz,
	double origx, double origy, double origz,
	double dx,  double dy,  double dz, double value, double integral, void *pdata);
OctNode* newOctNode(OctNode* parent, Vect3D center, Vect3D ds, double value, void *pdata);
void subdivide(OctNode* parent, double value, void *pdata);
void freeTree3D(Tree3D* tree, int freeDataPointer);
void freeNode(OctNode* node, int freeDataPointer);
OctNode* getLeafNodeContaining(Tree3D* tree, Vect3D coords);
OctNode* getLeafContaining(OctNode* node, double x, double y, double z);

ResultTreeNode* addResult(ResultTreeNode* prtn, double value, double volume, OctNode* pnode);
void freeResultTree(ResultTreeNode* prtn);
ResultTreeNode*  getHighestValue(ResultTreeNode* prtn);
ResultTreeNode* getHighestLeafValue(ResultTreeNode* prtree);
ResultTreeNode* getHighestLeafValueMinSize(ResultTreeNode* prtree, double sizeMinX, double sizeMinY, double sizeMinZ);

Tree3D* readTree3D(FILE *fpio);
int readNode(FILE *fpio, OctNode* node);
int writeTree3D(FILE *fpio, Tree3D* tree);
int writeNode(FILE *fpio, OctNode* node);

int nodeContains(OctNode* node, double x, double y, double z);
int extendedNodeContains(OctNode* node, double x, double y, double z, int checkZ);


/* */
/*------------------------------------------------------------/ */


