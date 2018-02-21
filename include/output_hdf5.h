//  Output hdf5
#include "hdf5.h"
//hid_t CreateH5_file_Para(char *FileName, MPI_Comm comm,MPI_Info info);
//void WriteData_Para(hid_t fileORgroupID, char *DatasetName, hid_t DataType, size_t Rank, size_t *BufDims, size_t *Offset, void *data);
hid_t CreateH5_file_Seq(char *FileName);
hid_t OpenH5_file(char *FileName);
hid_t CreateH5_group(hid_t fileID, char *GroupName);
void CreateH5_dataset(hid_t fileORgroupID, char *DatasetName, hid_t DataType , size_t Rank, size_t *DimInfo);
void WriteData_Seq(hid_t fileORgroupID, char *DatasetName,hid_t DataType, size_t Rank, size_t *BufDims, size_t *Offset, void *data);
