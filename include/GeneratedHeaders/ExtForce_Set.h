//This file is generated by Bash-script

int GAPS_APT_SetExtForceInfo(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	pPtc->num_forces=0;
	//Count force loaded
	if(0!=pInputs->ExtForce_Cal_RadLarmor)
	{
		(pPtc->num_forces)++;
	}
	if(0!=pInputs->ExtForce_Cal_GCElecCollision)
	{
		(pPtc->num_forces)++;
	}
	if(0!=pInputs->ExtForce_Cal_GCElecBremsstrahlung)
	{
		(pPtc->num_forces)++;
	}
	pPtc->ForceName= (char **)calloc((pPtc->num_forces),sizeof(char *));
	pPtc->ForceType= (int *)calloc((pPtc->num_forces),sizeof(int));

	//Set Force name and type
	int index=0;
	if(0!=pInputs->ExtForce_Cal_RadLarmor)
	{
		pPtc->ForceName[index] = "RadLarmor";
		pPtc->ForceType[index] = 0;
		index++;
	}
	if(0!=pInputs->ExtForce_Cal_GCElecCollision)
	{
		pPtc->ForceName[index] = "GCElecCollision";
		pPtc->ForceType[index] = 1;
		index++;
	}
	if(0!=pInputs->ExtForce_Cal_GCElecBremsstrahlung)
	{
		pPtc->ForceName[index] = "GCElecBremsstrahlung";
		pPtc->ForceType[index] = 2;
		index++;
	}
	return 0;
}