
#include "GeneratedMeshDeformersLibrary.h"
#include "DynamicMesh3.h"
#include "FrameTypes.h"
#include "MeshNormals.h"
#include "Async/ParallelFor.h"


UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshAxisSinWave1D(UGeneratedMesh* MeshObj, float Magnitude, float Frequency, float FrequencyShift, FVector AxisIn, FVector UpIn)
{
	FVector3d Axis(AxisIn), UpVector(UpIn);
	Axis.Normalize();
	UpVector.Normalize();

	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);
					double Dot = Pos.Dot(Axis);
					Dot += FrequencyShift;
					FVector3d NewPos = Pos + Magnitude * FMathd::Sin(Frequency * Dot) * UpVector;
					Mesh.SetVertex(vid, NewPos);
				}
			});
		});
	}

	return MeshObj;
}



UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshMove(UGeneratedMesh* MeshObj, FVector FromWorldLocation /*= FVector(0, 0, 0)*/, FVector ToWorldLocation /*= FVector(0, 0, 0)*/, float Radius /*= 1*/, float Hardness /*= 1*/, float Magnitude /*= 1*/)
{
	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);
					float Distance = (Pos - FromWorldLocation).Length();
					
					if (Distance < Radius)
					{
						float Ratio = FMath::Clamp((Radius - Distance) / Radius, Hardness, 1.f);	//hardness makes a bottom clamp
						float Scale = Ratio * Magnitude;	//magnitude scale whole operation
						//float Displacement = (ToWorldLocation - FromWorldLocation).Size();
						FVector3d NewPos = Pos + ((ToWorldLocation - FromWorldLocation) * Scale);
						Mesh.SetVertex(vid, NewPos);
					}
					
				}
			});
		});
	}
	return MeshObj;
}


float UGeneratedMeshDeformersLibrary::HeightAtPixel(float X, float Y, void* TexturePointer, int32 TextureHeight, int32 TextureWidth, int32 BytesPerPixel /*= 4*/)
{
	uint8* MipData = static_cast<uint8*>(TexturePointer);

	//We assume texture pointer is locked for read...
	int32 YRounded = FMath::RoundToInt(Y);
	int32 XRounded = FMath::RoundToInt(X);

	bool bValidSample = XRounded < TextureHeight && YRounded < TextureWidth;

	if (!bValidSample)
	{
		return 0;
	}

	//BytesPerPixel
	int32 SampleIndex = ((TextureWidth * YRounded) + XRounded) * BytesPerPixel;
	
	if (SampleIndex >= (BytesPerPixel * TextureHeight * TextureWidth))
	{
		return 0;
	}
	return MipData[SampleIndex] / 255.f;

	//Temp: get closest, TODO: interpolate from 4 nearest points
}

UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshHeightmap(UGeneratedMesh* MeshObj, UTexture2D* HeightmapTexture, float ZScale/*= 1*/, FVector DirectionNormal/*= FVector(0, 0, 1)*/)
{
	//Lock texture for read and get some metadata
	int32 BytesPerPixel = 4;

	EPixelFormat Format = HeightmapTexture->GetPixelFormat();
	if (Format != PF_B8G8R8A8)
	{
		UE_LOG(LogTemp, Log, TEXT("NB: Pixel format is not PF_R8G8B8A8: %d"), Format);
		//Todo: detect format and adjust BytesPerPixel

		BytesPerPixel = 1;
	}

	const int32 TextureWidth = HeightmapTexture->PlatformData->Mips[0].SizeX;
	const int32 TextureHeight= HeightmapTexture->PlatformData->Mips[0].SizeY;
	const int32 DataLength = TextureHeight * TextureWidth * BytesPerPixel;

	void* TextureDataPointer = HeightmapTexture->PlatformData->Mips[0].BulkData.Lock(LOCK_READ_ONLY);

	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			//grab bounds for mapped sampling
			//TODO: generalize to projection vector
			FAxisAlignedBox3d Bounds = Mesh.GetBounds();
			float MeshWidth = Bounds.Width();
			float MeshHeight = Bounds.Height();
			FVector3d Center = Bounds.Center();

			//UE_LOG(LogTemp, Log, TEXT("MeshSize: %1.3f, %1.3f"), MeshWidth, MeshHeight);

			//TODO: optimize loop
			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);

					//Normalize sample
					float XSample = (Pos.X + (MeshHeight / 2)) / MeshHeight * TextureHeight;
					float YSample = (Pos.Y + (MeshWidth / 2)) / MeshWidth * TextureWidth;

					//UE_LOG(LogTemp, Log, TEXT("Sample: %1.3f, %1.3f"), XSample, YSample);

					float UnscaledZ = HeightAtPixel(XSample, YSample, TextureDataPointer, TextureHeight, TextureWidth, BytesPerPixel);

					//TODO: support vector deform

					//Apply new delta
					FVector3d NewPos = Pos;
					NewPos.Z = Pos.Z + UnscaledZ * ZScale;

					Mesh.SetVertex(vid, NewPos);
				}
			});
		});
	}

	//re-lock heightmap
	HeightmapTexture->PlatformData->Mips[0].BulkData.Unlock();

	return MeshObj;
}

UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshAxisSinWaveRadial(UGeneratedMesh* MeshObj, float Magnitude, float Frequency, float FrequencyShift, FVector AxisIn)
{
	FVector3d Axis(AxisIn);
	Axis.Normalize();
	FFrame3d AxisFrame(FVector3d::Zero(), Axis);

	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);
					double Dot = Pos.Dot(Axis);
					Dot += FrequencyShift;
					double Displacement = Magnitude * FMathd::Sin(Frequency * (Dot + FrequencyShift));
					FVector3d PlaneVec = AxisFrame.ToPlane(Pos);
					PlaneVec.Normalize();
					FVector3d NewPos = Pos + Displacement * PlaneVec;
					Mesh.SetVertex(vid, NewPos);
				}
			});
		});
	}

	return MeshObj;
}



UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshPerlinNoiseNormal(UGeneratedMesh* MeshObj, float Magnitude, float Frequency, FVector FrequencyShift, int32 RandomSeed)
{
	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			FMath::SRandInit(RandomSeed);
			const float RandomOffset = 10000.0f * FMath::SRand();
			FVector3d Offset(RandomOffset, RandomOffset, RandomOffset);
			Offset += (FVector3d)FrequencyShift;

			FMeshNormals Normals(&Mesh);
			Normals.ComputeVertexNormals();

			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);
					FVector NoisePos = (FVector)( (double)Frequency * (Pos + Offset) );
					float Displacement = Magnitude * FMath::PerlinNoise3D(Frequency * NoisePos);
					FVector3d NewPos = Pos + Displacement * Normals[vid];
					Mesh.SetVertex(vid, NewPos);
				}
			});
		});
	}

	return MeshObj;
}




UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshPerlinAlongAxisNormal(UGeneratedMesh* MeshObj, float Magnitude /*= 1*/, float Frequency /*= 1*/, FVector FrequencyShift /*= FVector(0, 0, 0)*/, int32 RandomSeed /*= 31337*/, FVector DirectionNormal /*= FVector(0, 0, 1)*/)
{
	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
			{
				FMath::SRandInit(RandomSeed);
				const float RandomOffset = 10000.0f * FMath::SRand();
				FVector3d Offset(RandomOffset, RandomOffset, RandomOffset);
				Offset += (FVector3d)FrequencyShift;

				FMeshNormals Normals(&Mesh);
				Normals.ComputeVertexNormals();

				ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
					{
						if (Mesh.IsVertex(vid))
						{
							FVector3d Pos = Mesh.GetVertex(vid);
							FVector NoisePos = (FVector)((double)Frequency * (Pos + Offset));
							float Displacement = Magnitude * FMath::PerlinNoise3D(Frequency * NoisePos);
							FVector VertexNormal = (FVector)Normals[vid];

							if (acosf(VertexNormal | DirectionNormal) < (PI/4))
							{
								FVector3d NewPos = Pos + Displacement * VertexNormal;
								Mesh.SetVertex(vid, NewPos);
							}							
						}
					});
			});
	}

	return MeshObj;
}


UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshPerlinNormalPastHeight(UGeneratedMesh* MeshObj, float Magnitude /*= 1*/, float Frequency /*= 1*/, FVector FrequencyShift /*= FVector(0, 0, 0)*/, int32 RandomSeed /*= 31337*/, FVector HeightVector /*= FVector(0, 0, 1)*/)
{
	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
			{
				FMath::SRandInit(RandomSeed);
				const float RandomOffset = 10000.0f * FMath::SRand();
				FVector3d Offset(RandomOffset, RandomOffset, RandomOffset);
				Offset += (FVector3d)FrequencyShift;

				FMeshNormals Normals(&Mesh);
				Normals.ComputeVertexNormals();

				ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
					{
						if (Mesh.IsVertex(vid))
						{
							FVector3d Pos = Mesh.GetVertex(vid);
							FVector NoisePos = (FVector)((double)Frequency * (Pos + Offset));
							float Displacement = Magnitude * FMath::PerlinNoise3D(Frequency * NoisePos);
							FVector VertexNormal = (FVector)Normals[vid];

							if (Pos.Z >= HeightVector.Z)
							{
								FVector3d NewPos = Pos + Displacement * FVector(0,0,1);
								Mesh.SetVertex(vid, NewPos);
							}
						}
					});
			});
	}

	return MeshObj;
}

UGeneratedMesh* UGeneratedMeshDeformersLibrary::SmoothMeshUniform(UGeneratedMesh* MeshObj, float Alpha, int32 Iterations)
{
	Alpha = FMathf::Clamp(Alpha, 0.0f, 1.0f);
	Iterations = FMath::Clamp(Iterations, 0, 100);

	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			int32 NumV = Mesh.MaxVertexID();
			TArray<FVector3d> SmoothPositions;
			SmoothPositions.SetNum(NumV);

			for (int32 k = 0; k < Iterations; ++k)
			{
				ParallelFor(NumV, [&](int32 vid)
				{
					if (Mesh.IsVertex(vid))
					{
						FVector3d Centroid;
						Mesh.GetVtxOneRingCentroid(vid, Centroid);
						SmoothPositions[vid] = FVector3d::Lerp(Mesh.GetVertex(vid), Centroid, Alpha);
					}
				});
				for (int32 vid = 0; vid < NumV; ++vid)
				{
					if (Mesh.IsVertex(vid))
					{
						Mesh.SetVertex(vid, SmoothPositions[vid]);
					}
				}
			}
		});
	}

	return MeshObj;
}

void InitializeBrushIndices(int MapSize, int Radius, TArray<TArray<int32>>& ErosionBrushIndices, TArray<TArray<float>>& ErosionBrushWeights) {
	//ErosionBrushIndices = new int[mapSize * mapSize][];
	//ErosionBrushWeights = new float[mapSize * mapSize][];
	ErosionBrushIndices.SetNum(MapSize* MapSize);
	ErosionBrushWeights.SetNum(MapSize* MapSize);

	int* XOffsets = new int[Radius * Radius * 4];
	int* YOffsets = new int[Radius * Radius * 4];
	float* Weights = new float[Radius * Radius * 4];
	float WeightSum = 0;
	int AddIndex = 0;

	for (int i = 0; i < ErosionBrushIndices.Num(); i++) {
		int CentreX = i % MapSize;
		int CentreY = i / MapSize;

		if (CentreY <= Radius || CentreY >= MapSize - Radius || CentreX <= Radius + 1 || CentreX >= MapSize - Radius) {
			WeightSum = 0;
			AddIndex = 0;
			for (int Y = -Radius; Y <= Radius; Y++) {
				for (int X = -Radius; X <= Radius; X++) {
					float SqrDst = X * X + Y * Y;
					if (SqrDst < Radius * Radius) {
						int CoordX = CentreX + X;
						int CoordY = CentreY + Y;

						if (CoordX >= 0 && CoordX < MapSize && CoordY >= 0 && CoordY < MapSize) {
							float Weight = 1 - FMath::Sqrt(SqrDst) / Radius;
							WeightSum += Weight;
							Weights[AddIndex] = Weight;
							XOffsets[AddIndex] = X;
							YOffsets[AddIndex] = Y;
							AddIndex++;
						}
					}
				}
			}
		}

		int NumEntries = AddIndex;
		ErosionBrushIndices[i].SetNum(NumEntries);
		ErosionBrushWeights[i].SetNum(NumEntries);

		for (int j = 0; j < NumEntries; j++) {
			ErosionBrushIndices[i][j] = (YOffsets[j] + CentreY) * MapSize + XOffsets[j] + CentreX;
			ErosionBrushWeights[i][j] = Weights[j] / WeightSum;
		}
	}

	delete XOffsets, YOffsets, Weights;
}


FHeightAndGradient UGeneratedMeshDeformersLibrary::CalculateHeightAndGradient(TArray<float>& Nodes, int32 MapSize, float PosX, float PosY)
{
	int32 CoordX = (int32)PosX;
	int32 CoordY = (int32)PosY;

	// Calculate droplet's offset inside the cell (0,0) = at NW node, (1,1) = at SE node
	float X = PosX - CoordX;
	float Y = PosY - CoordY;

	// Calculate heights of the four nodes of the droplet's cell
	int32 NodeIndexNW = CoordY * MapSize + CoordX;
	float HeightNW = Nodes[NodeIndexNW];
	float HeightNE = Nodes[NodeIndexNW + 1];
	float HeightSW = Nodes[NodeIndexNW + MapSize];
	float HeightSE = Nodes[NodeIndexNW + MapSize + 1];

	// Calculate droplet's direction of flow with bilinear interpolation of height difference along the edges
	float GradientX = (HeightNE - HeightNW) * (1 - Y) + (HeightSE - HeightSW) * Y;
	float GradientY = (HeightSW - HeightNW) * (1 - X) + (HeightSE - HeightNE) * X;

	// Calculate height with bilinear interpolation of the heights of the nodes of the cell
	float Height = HeightNW * (1 - X) * (1 - Y) + HeightNE * X * (1 - Y) + HeightSW * (1 - X) * Y + HeightSE * X * Y;

	FHeightAndGradient Result;
	Result.Height = Height;
	Result.GradientX = GradientX;
	Result.GradientY = GradientY;

	return Result;
}

//Largely inspired from https://github.com/SebLague/Hydraulic-Erosion/blob/master/Assets/Scripts/Erosion.cs
UTexture2D* UGeneratedMeshDeformersLibrary::ErodeHeightMapTexture(UTexture2D* InTexture, const FHydroErosionParams& Params /*= FHydroErosionParams()*/)
{
	//only square textures
	if (!InTexture)
	{
		UE_LOG(LogTemp, Warning, TEXT("ErodeHeightMapTexture: Null texture passed for out or in, skipped."));
		return nullptr;
	}

	int32 SizeX = InTexture->PlatformData->Mips[0].SizeX;
	int32 SizeY = InTexture->PlatformData->Mips[0].SizeY;

	if (SizeX != SizeY)
	{
		UE_LOG(LogTemp, Warning, TEXT("ErodeHeightMapTexture: SizeX != SizeY, pass in square textures only, skipped."));
		return nullptr;
	}
	
	//Temp fixed params, todo: move to struct where we can tweak these params
	/*int32 Seed = 1;
	int32 ErosionRadius = 3;
	float Inertia = 0.05f; // At zero, water will instantly change direction to flow downhill. At 1, water will never change direction. 
	float SedimentCapacityFactor = 4; // Multiplier for how much sediment a droplet can carry
	float MinSedimentCapacity = 0.01f; // Used to prevent carry capacity getting too close to zero on flatter terrain
	float ErodeSpeed = 0.3f;
	float DepositSpeed = 0.3f;

	float EvaporateSpeed = 0.01f;
	float Gravity = 4;
	int32 MaxDropletLifetime = 30;*/

	float InitialWaterVolume = 1;
	float InitialSpeed = 1;

	int32 CurrentSeed;
	int32 CurrentErosionRadius;
	int32 CurrentMapSize;

	FRandomStream Prng = FRandomStream(Params.Seed);
	CurrentSeed = Params.Seed;
	CurrentErosionRadius = Params.ErosionRadius;

	TArray<TArray<int32>> ErosionBrushIndices;
	TArray<TArray<float>> ErosionBrushWeights;

	//Assuming it's one dim size
	int32 MapSize = InTexture->PlatformData->Mips[0].SizeX;

	InitializeBrushIndices(MapSize, Params.ErosionRadius, ErosionBrushIndices, ErosionBrushWeights);

	CurrentMapSize = MapSize;

	//TODO: get pixels to float pointer

	//inefficient: copy over into float array for first functioning implementation
	
	TArray<float> Map = Conv_GreyScaleTexture2DToFloatArray(InTexture);

	if (Map.Num() == 0)
	{
		UE_LOG(LogTemp, Warning, TEXT("ErodeHeightMapTexture Couldn't import texture, skipped."));
		return nullptr;
	}

	/*void* TextureDataPointer = InTexture->PlatformData->Mips[0].BulkData.Lock(LOCK_READ_ONLY);
	TArray<float> Map;
	Map.Reserve(MapSize * MapSize);

	for (int32 i = 0; i < MapSize; i++)
	{
		for (int32 j = 0; j < MapSize; j++)
		{
			float Height = HeightAtPixel(j, i, TextureDataPointer, MapSize, MapSize);
			Map.Add(Height);
		}
	}
	InTexture->PlatformData->Mips[0].BulkData.Unlock();*/
		
	//Start main algorithm loop
	for (int32 Iteration = 0; Iteration < Params.Iterations; Iteration++)
	{
		// Create water droplet at random point on map
		float PosX = Prng.FRandRange(0, MapSize - 1);
		float PosY = Prng.FRandRange(0, MapSize - 1);
		float DirX = 0;
		float DirY = 0;
		float Speed = InitialSpeed;
		float Water = InitialWaterVolume;
		float Sediment = 0;

		if (Iteration !=0 && Params.Iterations % Iteration == 5)
		{
			UE_LOG(LogTemp, Log, TEXT("Iteration: %d/%d"), Iteration, Params.Iterations);
		}

		for (int32 Lifetime = 0; Lifetime < Params.MaxDropletLifetime; Lifetime++)
		{
			int32 NodeX = (int)PosX;
			int32 NodeY = (int)PosY;
			int32 DropletIndex = NodeY * MapSize + NodeX;
			// Calculate droplet's offset inside the cell (0,0) = at NW node, (1,1) = at SE node
			float CellOffsetX = PosX - NodeX;
			float CellOffsetY = PosY - NodeY;

			// Calculate droplet's height and direction of flow with bilinear interpolation of surrounding heights
			FHeightAndGradient HeightAndGradient = CalculateHeightAndGradient(Map, MapSize, PosX, PosY);

			// Update the droplet's direction and position (move position 1 unit regardless of speed)
			DirX = (DirX * Params.Inertia - HeightAndGradient.GradientX * (1 - Params.Inertia));
			DirY = (DirY * Params.Inertia - HeightAndGradient.GradientY * (1 - Params.Inertia));
			// Normalize direction
			float Len = FMath::Sqrt(DirX * DirX + DirY * DirY);
			if (Len != 0) 
			{
				DirX /= Len;
				DirY /= Len;
			}
			PosX += DirX;
			PosY += DirY;

			// Stop simulating droplet if it's not moving or has flowed over edge of map
			if ((DirX == 0 && DirY == 0) || PosX < 0 || PosX >= MapSize - 1 || PosY < 0 || PosY >= MapSize - 1) 
			{
				break;
			}

			// Find the droplet's new height and calculate the deltaHeight
			float NewHeight = CalculateHeightAndGradient(Map, MapSize, PosX, PosY).Height;
			float DeltaHeight = NewHeight - HeightAndGradient.Height;

			// Calculate the droplet's sediment capacity (higher when moving fast down a slope and contains lots of water)
			float SedimentCapacity = FMath::Max(-DeltaHeight * Speed * Water * Params.SedimentCapacityFactor, Params.MinSedimentCapacity);

			// If carrying more sediment than capacity, or if flowing uphill:
			if (Sediment > SedimentCapacity || DeltaHeight > 0) 
			{
				// If moving uphill (deltaHeight > 0) try fill up to the current height, otherwise deposit a fraction of the excess sediment
				float AmountToDeposit = (DeltaHeight > 0) ? FMath::Min(DeltaHeight, Sediment) : (Sediment - SedimentCapacity) * Params.DepositSpeed;
				Sediment -= AmountToDeposit;

				// Add the sediment to the four nodes of the current cell using bilinear interpolation
				// Deposition is not distributed over a radius (like erosion) so that it can fill small pits
				Map[DropletIndex] += AmountToDeposit * (1 - CellOffsetX) * (1 - CellOffsetY);
				Map[DropletIndex + 1] += AmountToDeposit * CellOffsetX * (1 - CellOffsetY);
				Map[DropletIndex + MapSize] += AmountToDeposit * (1 - CellOffsetX) * CellOffsetY;
				Map[DropletIndex + MapSize + 1] += AmountToDeposit * CellOffsetX * CellOffsetY;

				if (Params.bPreviewDropletDepositionPaths)
				{
					Map[DropletIndex] = 1.0f;
					Map[DropletIndex + 1] = 0.9f;
					Map[DropletIndex + MapSize] = 0.9f;
					Map[DropletIndex + MapSize + 1] = 0.9f;
				}

				//UE_LOG(LogTemp, Log, TEXT("Output: #%d Deposit: %1.3f"), Iteration, AmountToDeposit);

			}
			else {
				// Erode a fraction of the droplet's current carry capacity.
				// Clamp the erosion to the change in height so that it doesn't dig a hole in the terrain behind the droplet
				float AmountToErode = FMath::Min((SedimentCapacity - Sediment) * Params.ErodeSpeed, -DeltaHeight);

				//UE_LOG(LogTemp, Log, TEXT("Output: #%d Erode: %1.3f"), Iteration, AmountToErode);

				// Use erosion brush to erode from all nodes inside the droplet's erosion radius
				for (int32 BrushPointIndex = 0; BrushPointIndex < ErosionBrushIndices[DropletIndex].Num(); BrushPointIndex++) 
				{
					int32 NodeIndex = ErosionBrushIndices[DropletIndex][BrushPointIndex];
					float WeighedErodeAmount = AmountToErode * ErosionBrushWeights[DropletIndex][BrushPointIndex];
					float DeltaSediment = (Map[NodeIndex] < WeighedErodeAmount) ? Map[NodeIndex] : WeighedErodeAmount;
					Map[NodeIndex] -= DeltaSediment;

					if (Params.bPreviewDropletErosionPaths)
					{
						Map[DropletIndex] = 0.f;
					}

					Sediment += DeltaSediment;
				}
			}

			// Update droplet's speed and water content
			Speed = FMath::Sqrt(Speed * Speed + DeltaHeight * Params.Gravity);
			Water *= (1 - Params.EvaporateSpeed);
		}
	}

	//Copy back data to out texture (todo: fast update/edit in place)
	return Conv_GrayScaleFloatArrayToTexture2D(Map, FVector2D(MapSize, MapSize));
}

UTexture2D* UGeneratedMeshDeformersLibrary::GenerateTransientCopy(UTexture2D* InTexture)
{
	if (!InTexture->IsValidLowLevelFast())
	{
		UE_LOG(LogTemp, Warning, TEXT("GenerateTransientCopy: Invalid Texture to copy."));
		return nullptr;
	}
	UTexture2D* TexturePointer = UTexture2D::CreateTransient(
		InTexture->PlatformData->Mips[0].SizeX, 
		InTexture->PlatformData->Mips[0].SizeY,
		InTexture->GetPixelFormat());
	TexturePointer->UpdateResource();
	return TexturePointer;
}

TArray<float> UGeneratedMeshDeformersLibrary::Conv_GreyScaleTexture2DToFloatArray(UTexture2D* InTexture)
{
	TArray<float> FloatArray;

	//sanity check for types we support atm
	if (InTexture->PlatformData->PixelFormat != PF_B8G8R8A8 &&
		InTexture->PlatformData->PixelFormat != PF_R8G8B8A8)
	{
		UE_LOG(LogTemp, Warning, TEXT("Invalid float array conversion not yet supported for requested pixel format."));
		return FloatArray;
	}

	FloatArray.SetNum(InTexture->GetSizeX() * InTexture->GetSizeY());

	// Lock the texture so it can be read
	uint8* MipData = static_cast<uint8*>(InTexture->PlatformData->Mips[0].BulkData.Lock(LOCK_READ_ONLY));

	for (int i = 0; i < FloatArray.Num(); i++)
	{
		int MipPointer = i * 4;
		float GreyscaleValue = (MipData[MipPointer] + MipData[MipPointer + 1] + MipData[MipPointer + 2]) / 3.f;
		FloatArray[i] = GreyscaleValue / 255.f;	 //normalize it
	}

	// Unlock the texture
	InTexture->PlatformData->Mips[0].BulkData.Unlock();

	return FloatArray;
}

UTexture2D* UGeneratedMeshDeformersLibrary::Conv_GrayScaleFloatArrayToTexture2D(const TArray<float>& InFloatArray, const FVector2D InSize /*= FVector2D(0, 0)*/)
{
	FVector2D Size;

	//Create square image and lock for writing
	if (InSize == FVector2D(0, 0))
	{
		int32 Length = FMath::Pow(InFloatArray.Num(), 0.5);
		if (Length * Length != InFloatArray.Num())
		{
			UE_LOG(LogTemp, Warning, TEXT("Invalid float array without specified size, needs to be square."));
			return nullptr;
		}
		Size = FVector2D(Length, Length);
	}
	else
	{
		Size = InSize;
	}

	UTexture2D* Pointer = UTexture2D::CreateTransient(Size.X, Size.Y, PF_B8G8R8A8);
	uint8* MipData = static_cast<uint8*>(Pointer->PlatformData->Mips[0].BulkData.Lock(LOCK_READ_WRITE));

	//Copy Data
	for (int i = 0; i < InFloatArray.Num(); i++)
	{
		int MipPointer = i * 4;
		int GreyValue = InFloatArray[i] * 255.f;
		MipData[MipPointer] = GreyValue;
		MipData[MipPointer + 1] = GreyValue;
		MipData[MipPointer + 2] = GreyValue;
		MipData[MipPointer + 3] = 255;	//Alpha
	}

	//Unlock and Return data
	Pointer->PlatformData->Mips[0].BulkData.Unlock();
	Pointer->UpdateResource();
	return Pointer;
}

