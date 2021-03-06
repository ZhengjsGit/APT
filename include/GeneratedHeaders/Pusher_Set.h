//This file is generated by Bash-script

int GAPS_APT_SetParticlePusher(Gaps_APT_ParticlePusher *pPusher,Gaps_IO_InputsContainer *pInputs)
{
	int Type = pInputs->Pusher_Type;
	if( Type<0 || Type >=13)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_SetParticlePusher: Pusher type is wrong! You input Push_Type = %ld\n",(long)Type);
	}
	else
	{
		if(0 == Type)
		{
			*pPusher = GAPS_APT_Pusher_RVPA_Cay3D;
		}
		if(1 == Type)
		{
			*pPusher = GAPS_APT_Pusher_LCCSA_SymEuler;
		}
		if(2 == Type)
		{
			*pPusher = GAPS_APT_Pusher_RCSA_SymEuler;
		}
		if(3 == Type)
		{
			*pPusher = GAPS_APT_Pusher_RVPA_Exp3D;
		}
		if(4 == Type)
		{
			*pPusher = GAPS_APT_Pusher_RungeKutta;
		}
		if(5 == Type)
		{
			*pPusher = GAPS_APT_Pusher_RNCSA_4D;
		}
		if(6 == Type)
		{
			*pPusher = GAPS_APT_Pusher_RECSA_GF4D;
		}
		if(7 == Type)
		{
			*pPusher = GAPS_APT_Pusher_LCCSA_IMP;
		}
		if(8 == Type)
		{
			*pPusher = GAPS_APT_Pusher_Regular_Boris;
		}
		if(9 == Type)
		{
			*pPusher = GAPS_APT_Pusher_Onehalf_Boris;
		}
		if(10 == Type)
		{
			*pPusher = GAPS_APT_Pusher_CSA_SymEuler;
		}
		if(11 == Type)
		{
			*pPusher = GAPS_APT_Pusher_CSA_imEuler;
		}
		if(12 == Type)
		{
			*pPusher = GAPS_APT_Pusher_Stormer_Verlet;
		}
	}
	return 0;
}
