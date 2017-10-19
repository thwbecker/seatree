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


/*   octree.c

	octree search functions

*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
	history:

	ver 01    02NOV2000  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#define EXTERN_MODE 1

#include "geometry.h"
#include "octtree.h"
#include "ran1.h"


/*** function to create a new OctNode */

OctNode* newOctNode(OctNode* parent, Vect3D center, Vect3D ds, double value, void *pdata)
{

	int l, m, n;
	OctNode* node;

	node = (OctNode*) malloc(sizeof(OctNode));

	node->parent = parent;
	node->center = center;
	node->ds = ds;
	node->value = value;
	node->pdata = pdata;

	for (l = 0; l < 2; l++)
		for (m = 0; m < 2; m++)
			for (n = 0; n < 2; n++)
		node->child[l][m][n] = NULL;

	node->isLeaf = 1;

	return(node);
}



/*** function to create a new Tree3D - an x, y, z array of octtree root nodes ***/

Tree3D* newTree3D(int data_code, int numx, int numy, int numz,
	double origx, double origy, double origz,
	double dx,  double dy,  double dz, double value, double integral, void *pdata)
{

	int ix, iy, iz;
	OctNode**** garray;
	Tree3D* tree;
	Vect3D center, ds;

	// alocate nodes in x, y and z

	if ((garray = (OctNode ****) malloc((size_t) numx * sizeof(OctNode***))) == NULL)
		return(NULL);

	ds.x = dx;
	ds.y = dy;
	ds.z = dz;

	for (ix = 0; ix < numx; ix++) {
		center.x = origx + (double) ix * dx + dx / 2.0;
        	if ((garray[ix] = (OctNode ***) malloc((size_t) numy *
				sizeof(OctNode**))) == NULL)
			return(NULL);
		for (iy = 0; iy < numy; iy++) {
			center.y = origy + (double) iy * dy + dy / 2.0;
       			if ((garray[ix][iy] = (OctNode **) malloc((size_t) numz *
					sizeof(OctNode*))) == NULL)
				return(NULL);
			for (iz = 0; iz < numz; iz++) {
				center.z = origz + (double) iz * dz + dz / 2.0;
				garray[ix][iy][iz] = newOctNode(NULL, center, ds, value, pdata);
			}
		}
	}


	// create Tree3D structure

	tree = (Tree3D*) malloc(sizeof(Tree3D));

	tree->nodeArray = garray;
	tree->data_code = data_code;
	tree->numx = numx;
	tree->numy = numy;
	tree->numz = numz;
	tree->orig.x = origx;
	tree->orig.y = origy;
	tree->orig.z = origz;
	tree->ds = ds;
	tree->integral = integral;

	return(tree);

}



/*** function to sudivide a node into child nodes ***/

void subdivide(OctNode* parent, double value, void *pdata) {

	int ix, iy, iz;
	Vect3D center, ds;

	ds.x = parent->ds.x / 2.0;
	ds.y = parent->ds.y / 2.0;
	ds.z = parent->ds.z / 2.0;

	// alocate nodes in x, y and z

	for (ix = 0; ix < 2; ix++) {
		center.x = parent->center.x + (double) (2 * ix - 1) * ds.x / 2.0;
		for (iy = 0; iy < 2; iy++) {
			center.y = parent->center.y + (double) (2 * iy - 1) * ds.y / 2.0;
			for (iz = 0; iz < 2; iz++) {
				center.z = parent->center.z +
					(double) (2 * iz - 1) * ds.z / 2.0;
				parent->child[ix][iy][iz] =
					newOctNode(parent, center, ds, value, pdata);
			}
		}
	}

	if (parent != NULL)
		parent->isLeaf = 0;

}


/*** function to free a Tree3D ***/

void freeTree3D(Tree3D* tree, int freeDataPointer)
{

	int ix, iy, iz;

	for (ix = 0; ix < tree->numx; ix++) {
		for (iy = 0; iy < tree->numy; iy++) {
			for (iz = 0; iz < tree->numz; iz++) {
				freeNode(tree->nodeArray[ix][iy][iz], freeDataPointer);
			}
       			free(tree->nodeArray[ix][iy]);
		}
        	free(tree->nodeArray[ix]);
	}

	free(tree);

}



/*** function to free an OctNode and all its child nodes ***/

void freeNode(OctNode* node, int freeDataPointer) {

	int ix, iy, iz;
	for (ix = 0; ix < 2; ix++) {
		for (iy = 0; iy < 2; iy++) {
			for (iz = 0; iz < 2; iz++) {
				if (node->child[ix][iy][iz] != NULL)
					freeNode(node->child[ix][iy][iz], freeDataPointer);
			}
		}
	}

	// try to free data
	if (freeDataPointer)
		free(node->pdata);
	free(node);

}


/*** function to get the leaf node in a Tree3D containing the given x, y, z coordinates ***/

OctNode* getLeafNodeContaining(Tree3D* tree, Vect3D coords)
{

	int ix, iy, iz;
	OctNode* rootNode;

	// get indices of root in tree
	ix = (int) ((coords.x - tree->orig.x) / tree->ds.x);
	if (ix < 0 || ix >= tree->numx)
		return(NULL);
	iy = (int) ((coords.y - tree->orig.y) / tree->ds.y);
	if (iy < 0 || iy >= tree->numy)
		return(NULL);
	iz = (int) ((coords.z - tree->orig.z) / tree->ds.z);
	if (iz < 0 || iz >= tree->numz)
		return(NULL);

	rootNode = tree->nodeArray[ix][iy][iz];

	if (rootNode == NULL)
		return(NULL);

	return(getLeafContaining(rootNode, coords.x, coords.y, coords.z));
}


/*** function to get the leaf node in a node containing the given x, y, z coordinates ***/

OctNode* getLeafContaining(OctNode* node, double x, double y, double z)
{

	int ix = 0, iy = 0, iz = 0;
	OctNode* childNode;

	// get indices of child in this node
	if (x >= node->center.x)
		ix = 1;
	if (y >= node->center.y)
		iy = 1;
	if (z >= node->center.z)
		iz = 1;

	childNode = node->child[ix][iy][iz];

	if (childNode == NULL)
		return(node);

	return(getLeafContaining(childNode, x, y, z));

}


/*** function to check if node contains the given x, y, z coordinates ***/

int nodeContains(OctNode* node, double x, double y, double z)
{

	Vect3D ds;
	double dx2, dy2, dz2;

	ds = node->ds;

	dx2 = ds.x / 2.0;
	if (x < node->center.x - dx2 || x > node->center.x + dx2)
		return(0);

	dy2 = ds.y / 2.0;
	if (y < node->center.y - dy2 || y > node->center.y + dy2)
		return(0);

	dz2 = ds.z / 2.0;
	if (z < node->center.z - dz2 || z > node->center.z + dz2)
		return(0);

	return(1);

}


/*** function to check if node or half-segments of adjacent equal-size nodes contains the given x, y, z coordinates ***/

int extendedNodeContains(OctNode* node, double x, double y, double z, int checkZ)
{

	Vect3D ds;
	double dx2, dy2, dz2;

	ds = node->ds;

	dx2 = ds.x;
	if (x < node->center.x - dx2 || x > node->center.x + dx2)
		return(0);

	dy2 = ds.y;
	if (y < node->center.y - dy2 || y > node->center.y + dy2)
		return(0);

	if (checkZ) {
		dz2 = ds.z;
		if (z < node->center.z - dz2 || z > node->center.z + dz2)
			return(0);
	}

	return(1);

}


/*** function to put Octtree node in results tree in order of value */

ResultTreeNode* addResult(ResultTreeNode* prtree, double value, double volume, OctNode* pnode)
{
	/* put address in result tree based on value */

	if (prtree == NULL) {	/* at empty node */
		if ((prtree = (ResultTreeNode* ) malloc(sizeof(ResultTreeNode))) == NULL)
			fprintf(stderr, "ERROR allocating memory for result-tree node.\n");
		prtree->value = value;
		prtree->volume = volume;	// node volume depends on geometry in physical space, may not be dx*dy*dz
		prtree->pnode = pnode;
		prtree->left = prtree->right = NULL;

	} else if (value == prtree->value)  {	// prevent assymetric tree if multiple identical values
		if (get_rand_int(-10000, 9999) < 0)
			prtree->left = addResult(prtree->left, value, volume, pnode);
		else
			prtree->right = addResult(prtree->right, value, volume, pnode);

	} else if (value < prtree->value)  {
		prtree->left = addResult(prtree->left, value, volume, pnode);

	} else  {
		prtree->right = addResult(prtree->right, value, volume, pnode);
	}

	return (prtree);
}


/*** function to free results tree */

void freeResultTree(ResultTreeNode* prtree)
{

	if (prtree->left != NULL)
		freeResultTree(prtree->left);
	if (prtree->right != NULL)
		freeResultTree(prtree->right);
	free(prtree);
}


/*** function to get ResultTreeNode with highest value */

ResultTreeNode* getHighestValue(ResultTreeNode* prtree)
{
	if (prtree->right == NULL)	/* right child is empty */
		return(prtree);

	return(getHighestValue(prtree->right));
}


/*** function to get ResultTree Leaf Node with highest value */

ResultTreeNode* getHighestLeafValue(ResultTreeNode* prtree)
{

	ResultTreeNode* prtree_returned = NULL;

	if (prtree->right != NULL)
		prtree_returned = getHighestLeafValue(prtree->right);

	if (prtree_returned != NULL)		// right leaf descendent
		return(prtree_returned);	// thus highest value leaf
	else					// no right leaf descendents
		if (prtree->pnode->isLeaf)	// this is leaf
			return(prtree);		// thus highest value leaf

	if (prtree->left != NULL)	// look for left leaf descendents
		prtree_returned = getHighestLeafValue(prtree->left);

	return(prtree_returned);
}



/*** function to get ResultTree Leaf Node with highest value */

ResultTreeNode* getHighestLeafValueMinSize(ResultTreeNode* prtree, double sizeMinX, double sizeMinY, double sizeMinZ)
{

	ResultTreeNode* prtree_returned = NULL;

	if (prtree->right != NULL)
		prtree_returned = getHighestLeafValueMinSize(prtree->right, sizeMinX, sizeMinY, sizeMinZ);

	if (prtree_returned != NULL)		// right leaf descendent
		return(prtree_returned);	// thus highest value leaf
	else {					// no right leaf descendents
		if (prtree->pnode->isLeaf	// this is leaf
			&& prtree->pnode->ds.x >= sizeMinX
			&& prtree->pnode->ds.y >= sizeMinY
			&& prtree->pnode->ds.z >= sizeMinZ) // not too small
				return(prtree);	// thus highest value leaf that is not too small
	}

	if (prtree->left != NULL)	// look for left leaf descendents
		prtree_returned = getHighestLeafValueMinSize(prtree->left, sizeMinX, sizeMinY, sizeMinZ);

	return(prtree_returned);
}






/*** function to read a Tree3D ***/

Tree3D* readTree3D(FILE *fpio)
{

	int istat;
	int istat_cum;
	int ix, iy, iz;

	Tree3D* tree = NULL;
	int data_code;
	int numx, numy, numz;
	Vect3D orig;
	Vect3D ds;
	double integral;


	istat = fread(&data_code, sizeof(int), 1, fpio);
	istat += fread(&numx, sizeof(int), 1, fpio);
	istat += fread(&numy, sizeof(int), 1, fpio);
	istat += fread(&numz, sizeof(int), 1, fpio);
	istat += fread(&orig, sizeof(Vect3D), 1, fpio);
	istat += fread(&ds, sizeof(Vect3D), 1, fpio);
	istat += fread(&integral, sizeof(double), 1, fpio);

	if (istat < 7)
		return(NULL);

	tree = newTree3D(data_code, numx, numy, numz, orig.x, orig.y, orig.z, ds.x, ds.y, ds.z, -1.0, integral, NULL);

	istat_cum = 0;
	for (ix = 0; ix < tree->numx; ix++) {
		for (iy = 0; iy < tree->numy; iy++) {
			for (iz = 0; iz < tree->numz; iz++) {
				istat += readNode(fpio, tree->nodeArray[ix][iy][iz]);
					if (istat < 0)
						return(NULL);
					istat_cum += istat;
			}
		}
	}

	return(tree);

}



/*** function to read an OctNode and all its child nodes ***/

int readNode(FILE *fpio, OctNode* node)
{

	int istat;
	int istat_cum;
	int ix, iy, iz;

	float value;


	istat = fread(&value, sizeof(float), 1, fpio);	/* node value */
	node->value = (double) value;
	istat += fread(&(node->isLeaf), sizeof(char), 1, fpio);		/* leaf flag, 1=leaf */

	if (istat < 2)
		return(-1);

	if (node->isLeaf)
		return(1);

	subdivide(node, -1.0, NULL);

	istat_cum = 1;

	for (ix = 0; ix < 2; ix++) {
		for (iy = 0; iy < 2; iy++) {
			for (iz = 0; iz < 2; iz++) {
				if (node->child[ix][iy][iz] != NULL) {
					istat = readNode(fpio, node->child[ix][iy][iz]);
					if (istat < 0)
						return(-1);
					istat_cum += istat;
				}
			}
		}
	}

	return(istat_cum);

}



/*** function to write a Tree3D ***/

int writeTree3D(FILE *fpio, Tree3D* tree)
{

	int istat;
	int istat_cum;
	int ix, iy, iz;

	istat = fwrite(&(tree->data_code), sizeof(int), 1, fpio);
	istat += fwrite(&(tree->numx), sizeof(int), 1, fpio);
	istat += fwrite(&(tree->numy), sizeof(int), 1, fpio);
	istat += fwrite(&(tree->numz), sizeof(int), 1, fpio);
	istat += fwrite(&(tree->orig), sizeof(Vect3D), 1, fpio);
	istat += fwrite(&(tree->ds), sizeof(Vect3D), 1, fpio);
	istat += fwrite(&(tree->integral), sizeof(double), 1, fpio);

	if (istat < 6)
		return(-1);

	istat_cum = 0;
	for (ix = 0; ix < tree->numx; ix++) {
		for (iy = 0; iy < tree->numy; iy++) {
			for (iz = 0; iz < tree->numz; iz++) {
				istat = writeNode(fpio, tree->nodeArray[ix][iy][iz]);
				if (istat < 0)
					return(-1);
				istat_cum += istat;
			}
		}
	}

	return(istat_cum);

}



/*** function to write an OctNode and all its child nodes ***/

int writeNode(FILE *fpio, OctNode* node)
{

	int istat;
	int istat_cum;
	int ix, iy, iz;

	float value;

	value = (float) node->value;
	istat = fwrite(&value, sizeof(float), 1, fpio);	/* node value */
	istat += fwrite(&(node->isLeaf), sizeof(char), 1, fpio);		/* leaf flag, 1=leaf */

	if (istat < 2)
		return(-1);

	if (node->isLeaf)
		return(1);

	istat_cum = 1;
	for (ix = 0; ix < 2; ix++) {
		for (iy = 0; iy < 2; iy++) {
			for (iz = 0; iz < 2; iz++) {
				if (node->child[ix][iy][iz] != NULL) {
					istat = writeNode(fpio, node->child[ix][iy][iz]);
					if (istat < 0)
						return(-1);
					istat_cum += istat;
				}
			}
		}
	}

	return(istat_cum);

}




/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

