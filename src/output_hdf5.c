#ifdef DATA_FORMAT_HDF5
#include "output_hdf5.h"

hid_t CreateH5_file_Seq(char *FileName)
{
	return H5Fcreate (FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

/*hid_t CreateH5_file_Para(char *FileName, MPI_Comm comm,MPI_Info info)
{
	hid_t fileID, plist_id;

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id,comm,info);

	fileID =  H5Fcreate (FileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	
	H5Pclose(plist_id);
	return fileID;
}
*/

hid_t OpenH5_file(char *FileName)
{
	return H5Fopen(FileName, H5F_ACC_RDWR, H5P_DEFAULT);
}

hid_t CreateH5_group(hid_t fileID, char *GroupName)
{
	hid_t groupID;

	groupID=H5Gcreate2(fileID,GroupName,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	return groupID;
}

void CreateH5_dataset(hid_t fileORgroupID, char *DatasetName, hid_t DataType , size_t Rank, size_t *DimInfo)
{//DimInfo is a 1D-array: DimInfo[rank+1]={rank, num_1, num_2, ... , num_rank}, rank is the dimention of data, num_1~num_rank is the number of each dimention.
	
	hid_t filespace, datasetID;
	hsize_t dimsf[Rank];

	int i;
	for(i=0;i<Rank;i++)
	{
		(*(dimsf+i)) = (*(DimInfo+i));
	}

	filespace = H5Screate_simple(Rank,dimsf,NULL);

	datasetID = H5Dcreate(fileORgroupID,DatasetName,DataType,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

	H5Sclose(filespace);
	H5Dclose(datasetID);
}

/*void WriteData_Para(hid_t fileORgroupID, char *DatasetName, hid_t DataType, size_t Rank, size_t *BufDims, size_t *Offset, void *data)
{
	herr_t status;
	hid_t datasetID = H5Dopen2(fileORgroupID, DatasetName, H5P_DEFAULT);	

	hsize_t count[Rank], offset[Rank];

	int i;
	for (i=0;i<Rank;i++)
	{
		(*(count+i)) = (*(BufDims+i));
		(*(offset+i)) = (*(Offset+i));
	}

	hid_t filespace, memspace, plist_id;
	memspace = H5Screate_simple(Rank,count,NULL);

	filespace = H5Dget_space(datasetID);
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET,offset,NULL,count,NULL);

	plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	status = H5Dwrite(datasetID, DataType,memspace,filespace,plist_id,data);

	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Pclose(plist_id);
	H5Dclose(datasetID);
}
*/

void WriteData_Seq(hid_t fileORgroupID, char *DatasetName,hid_t DataType, size_t Rank, size_t *BufDims, size_t *Offset, void *data)
{
	herr_t status;
	hid_t datasetID = H5Dopen2(fileORgroupID, DatasetName, H5P_DEFAULT);	
	hsize_t count[Rank], offset[Rank];

	int i;
	for (i=0;i<Rank;i++)
	{
		(*(count+i)) = (*(BufDims+i));
		(*(offset+i)) = (*(Offset+i));
	}

	hid_t filespace, memspace;
	memspace = H5Screate_simple(Rank,count,NULL);

	filespace = H5Dget_space(datasetID);
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET,offset,NULL,count,NULL);


	status = H5Dwrite(datasetID, DataType,memspace,filespace,H5P_DEFAULT,data);

	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Dclose(datasetID);
}

#endif
