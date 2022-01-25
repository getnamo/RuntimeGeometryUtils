#include "HeightmapDeformersLibrary.h"

float UHeightmapDeformersLibrary::HeightAtPixel(float X, float Y, void* TexturePointer, int32 TextureHeight, int32 TextureWidth, int32 BytesPerPixel /*= 4*/)
{
	uint8* MipData = static_cast<uint8*>(TexturePointer);

	//We assume texture pointer is locked for read...
	int32 YRounded = FMath::RoundToInt(Y);
	int32 XRounded = FMath::RoundToInt(X);

	bool bValidSample = XRounded < TextureHeight&& YRounded < TextureWidth;

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


float UHeightmapDeformersLibrary::HeightInSquareArray(const FTransform& SamplingTransform, TArray<float>& Array)
{
	int32 YRounded = FMath::RoundToInt(SamplingTransform.GetLocation().Y);
	int32 XRounded = FMath::RoundToInt(SamplingTransform.GetLocation().X);

	//TODO modify X/Y by transform rotation and scale

	int32 MapSize = (int32)FMath::Sqrt((float)Array.Num());

	bool bValidSample = XRounded < MapSize && YRounded < MapSize;

	if (!bValidSample)
	{
		UE_LOG(LogTemp, Warning, TEXT("UHeightmapDeformersLibrary::HeightInSquareArray Invalid sample 0.f returned."));
		return 0.f;
	}

	int32 SampleIndex = ((MapSize * YRounded) + XRounded);

	return Array[SampleIndex];
}


void UHeightmapDeformersLibrary::HydraulicErosionOnHeightMapWithInterrupt(TArray<float>& Map, const FHydroErosionParams& Params, TFunction<bool()>InterruptHandler /*= nullptr*/)
{
	if (Map.Num() == 0)
	{
		UE_LOG(LogTemp, Warning, TEXT("UHeightmapDeformersLibrary::HydraulicErosionOnHeightMap Empty map array passed."));
		return;
	}
	float InitialWaterVolume = 1;
	float InitialSpeed = 1;

	FRandomStream Prng = FRandomStream(Params.Seed);

	TArray<TArray<int32>> ErosionBrushIndices;
	TArray<TArray<float>> ErosionBrushWeights;

	//this should always be a square texture

	//todo: detect non-square texture by lack of square result
	int32 MapSize = (int32)FMath::Sqrt((float)Map.Num());

	//Todo: re-use this if possible
	InitializeBrushIndices(MapSize, Params.ErosionRadius, ErosionBrushIndices, ErosionBrushWeights);

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

		if (Iteration != 0 && Iteration % 10000 == 0)
		{
			UE_LOG(LogTemp, Log, TEXT("Iteration: %d/%d"), Iteration, Params.Iterations);
			if (InterruptHandler != nullptr)
			{
				if (InterruptHandler())
				{
					UE_LOG(LogTemp, Warning, TEXT("UHeightmapDeformersLibrary::HydraulicErosionOnHeightMapWithInterrupt Interrupted by handler"));
					return;
				}
			}
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
}

//From: https://github.com/svaarala/duktape/blob/master/misc/splitmix64.c
int64 UHeightmapDeformersLibrary::SplitMix64(int64& Seed)
{
	uint64_t x = Seed;
	uint64_t z = (x += UINT64_C(0x9E3779B97F4A7C15));
	z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
	z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
	Seed = x;
	return z ^ (z >> 31);
}

void UHeightmapDeformersLibrary::PerlinDeformMap(UPARAM(ref) TArray<float>& Map,
	float Magnitude /*= 1*/,
	float Frequency /*= 1*/, 
	FVector FrequencyShift /*= FVector(0, 0, 0)*/, 
	int32 RandomSeed /*= 31337*/, 
	int32 Octaves /*= 1*/,
	float OctaveFactor /*= 2.f*/,
	bool bRidged /*=false*/)
{
	FVector2D NoisePos;

	int32 MapSize = (int32)FMath::Sqrt((float)Map.Num());

	//Set X/Y from loop over square texture
	float OctaveFrequency = Frequency;
	float OctaveMagnitude = Magnitude;
	float Displacement = 0.f;

	for (int32 i = 0; i < Octaves; i++)
	{
		for (float Y = 0.f; Y < MapSize; Y++)
		{
			for (float X = 0.f; X < MapSize; X++)
			{
				int32 Index = FMath::RoundToInt((Y * MapSize) + X);

				NoisePos.X = X + FrequencyShift.X;
				NoisePos.Y = Y + FrequencyShift.Y;
				

				if (bRidged)
				{
					Displacement = OctaveMagnitude / 2.f * (1.f-FMath::Abs(FMath::PerlinNoise2D(NoisePos * OctaveFrequency)));
				}
				else
				{
					Displacement = OctaveMagnitude * (FMath::PerlinNoise2D(NoisePos * OctaveFrequency) + 1.f) / 2.f;
				}

				if (Index < Map.Num())
				{
					Map[Index] += Displacement;
				}
				else
				{
					UE_LOG(LogTemp, Warning, TEXT("Invalid index: %d/%d"), Index, Map.Num());
				}
			}
		}

		//Apply octave shift
		OctaveFrequency = OctaveFrequency * OctaveFactor;
		OctaveMagnitude = OctaveMagnitude / OctaveFactor;
	}
}


void UHeightmapDeformersLibrary::PerlinDeformMeshAlongCenter(
	UPARAM(ref) TArray<FVector>& InOutVertices,
	FVector Center,
	float Magnitude /*= 1*/,
	float Frequency /*= 1*/,
	FVector FrequencyShift /*= FVector(0, 0, 0)*/,
	int32 RandomSeed /*= 31337*/,
	int32 Octaves /*= 1*/,
	float OctaveFactor /*= 2.f*/,
	bool bRidged /*= false*/)
{
	FVector NoisePos;

	//Set X/Y from loop over square texture
	float OctaveFrequency = Frequency;
	float OctaveMagnitude = Magnitude;
	float Displacement = 0.f;

	//obtain average center for normal
	/*FVector Center = FVector(0.f);
	for (FVector& Vertex : InOutVertices)
	{
		Center += Vertex;
	}
	Center = Center / InOutVertices.Num();*/

	//shift every vertex by displacement
	for (int32 i = 0; i < Octaves; i++)
	{
		for (FVector& Vertex : InOutVertices)
		{
			NoisePos = Vertex + FrequencyShift;

			if (bRidged)
			{
				Displacement = OctaveMagnitude / 2.f * (1.f - FMath::Abs(FMath::PerlinNoise3D(NoisePos * OctaveFrequency)));
			}
			else
			{
				Displacement = OctaveMagnitude * (FMath::PerlinNoise3D(NoisePos * OctaveFrequency) + 1.f) / 2.f;
			}

			FVector Normal = (Vertex - Center).GetSafeNormal();

			Vertex = Vertex + (Displacement * Normal);	//should be along normal...
		}

		//Apply octave shift
		OctaveFrequency = OctaveFrequency * OctaveFactor;
		OctaveMagnitude = OctaveMagnitude / OctaveFactor;
	}
}

TArray<float> UHeightmapDeformersLibrary::SquareFloatMapSized(int32 OneSideLength)
{
	TArray<float>Map;
	Map.SetNumZeroed(OneSideLength * OneSideLength);

	return Map;
}


UTexture2D* UHeightmapDeformersLibrary::SquareTextureSized(int32 OneSideLength, EPixelFormat Format)
{
	UTexture2D* Pointer = UTexture2D::CreateTransient(OneSideLength, OneSideLength, Format);
	if (Format == EPixelFormat::PF_FloatRGBA)
	{
		Pointer->SRGB = false;
	}
	Pointer->UpdateResource();
	return Pointer;
}


void UHeightmapDeformersLibrary::CopyFloatArrayToTexture(const TArray<float>& SrcData, UTexture2D* TargetTexture)
{
	uint8* MipData = static_cast<uint8*>(TargetTexture->PlatformData->Mips[0].BulkData.Lock(LOCK_READ_WRITE));

	if (TargetTexture->GetPixelFormat() == EPixelFormat::PF_B8G8R8A8)
	{
		//Copy Data
		for (int i = 0; i < SrcData.Num(); i++)
		{
			int MipPointer = i * 4;
			int GreyValue = FMath::Clamp(SrcData[i], 0.f, 1.f) * 255.f;
			MipData[MipPointer] = GreyValue;
			MipData[MipPointer + 1] = GreyValue;
			MipData[MipPointer + 2] = GreyValue;
			MipData[MipPointer + 3] = 255;	//Alpha
		}
	}
	else if (TargetTexture->GetPixelFormat() == EPixelFormat::PF_G16)
	{
		for (int i = 0; i < SrcData.Num(); i++)
		{
			int32 MipPointer = i * 2;
			uint16 GreyValue = FMath::RoundToInt(FMath::Clamp(SrcData[i], 0.f, 1.f) * 65536.f);

			FMemory::Memcpy(&MipData[MipPointer], &GreyValue, sizeof(uint16));
		}
	}
	else if (TargetTexture->GetPixelFormat() == EPixelFormat::PF_FloatRGBA)
	{
		for (int i = 0; i < SrcData.Num(); i++)
		{
			int32 MipPointer = i * 8;
			FFloat16Color Color;
			
			Color.R = FMath::Clamp(SrcData[i], 0.f, 1.f);
			Color.G = Color.R;
			Color.B = Color.R;
			Color.A = 1.f;

			FMemory::Memcpy(&MipData[MipPointer], &Color, sizeof(FFloat16Color));
		}
	}

	//Unlock and Return data
	TargetTexture->PlatformData->Mips[0].BulkData.Unlock();
	TargetTexture->UpdateResource();
}

void UHeightmapDeformersLibrary::InitializeBrushIndices(int MapSize, int Radius, TArray<TArray<int32>>& ErosionBrushIndices, TArray<TArray<float>>& ErosionBrushWeights) {
	//ErosionBrushIndices = new int[mapSize * mapSize][];
	//ErosionBrushWeights = new float[mapSize * mapSize][];
	ErosionBrushIndices.SetNum(MapSize * MapSize);
	ErosionBrushWeights.SetNum(MapSize * MapSize);

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


FHeightAndGradient UHeightmapDeformersLibrary::CalculateHeightAndGradient(TArray<float>& Nodes, int32 MapSize, float PosX, float PosY)
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
UTexture2D* UHeightmapDeformersLibrary::HydraulicErosionOnHeightTexture(UTexture2D* InTexture, const FHydroErosionParams& Params /*= FHydroErosionParams()*/)
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

	//inefficient: copy over into float array for a first functioning implementation
	TArray<float> Map = Conv_GreyScaleTexture2DToFloatArray(InTexture);

	if (Map.Num() == 0)
	{
		UE_LOG(LogTemp, Warning, TEXT("ErodeHeightMapTexture Couldn't import texture, skipped."));
		return nullptr;
	}

	//Run the erosion algorithm
	HydraulicErosionOnHeightMap(Map, Params);

	//Copy back data to out texture (todo: fast update/edit in place)
	return Conv_GrayScaleFloatArrayToTexture2D(Map, FVector2D(SizeX, SizeY));
}


void UHeightmapDeformersLibrary::HydraulicErosionOnHeightMap(UPARAM(ref) TArray<float>& Map, const FHydroErosionParams& Params)
{
	HydraulicErosionOnHeightMapWithInterrupt(Map, Params, nullptr);
}

UTexture2D* UHeightmapDeformersLibrary::GenerateTransientCopy(UTexture2D* InTexture)
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


void UHeightmapDeformersLibrary::VectorOpArray(UPARAM(ref) TArray<float>& InOutArray, const TArray<float>& Other, EFloatAppendTypes AppendType /*= EFloatAppendTypes::Add*/, bool bNormalize /*= true*/)
{
	if (InOutArray.Num() != Other.Num())
	{
		UE_LOG(LogTemp, Warning, TEXT("UHeightmapDeformersLibrary::AppendAndRescale InOutArray and other size don't match. Skipped."))
			return;
	}

	float Max = 0;
	float Min = FLT_MAX;

	//Loop one apply change
	for (int32 i = 0; i < InOutArray.Num(); i++)
	{
		if (bNormalize)
		{
			if (InOutArray[i] < Min)
			{
				Min = InOutArray[i];
			}
			if (InOutArray[i] > Max)
			{
				Max = InOutArray[i];
			}
		}
		if (AppendType == EFloatAppendTypes::Add)
		{
			InOutArray[i] += Other[i];
		}
		else if (AppendType == EFloatAppendTypes::Subtract)
		{
			InOutArray[i] -= Other[i];
		}
		else if (AppendType == EFloatAppendTypes::Multiply)
		{
			InOutArray[i] *= Other[i];
		}
		else if (AppendType == EFloatAppendTypes::Divide)
		{
			InOutArray[i] /= Other[i];
		}
		else if (AppendType == EFloatAppendTypes::Override)
		{
			InOutArray[i] = Other[i];
		}
	}

	//loop 2, normalize result
	if (bNormalize)
	{
		float Range = Max - Min;

		for (int32 i = 0; i < InOutArray.Num(); i++)
		{
			InOutArray[i] = (InOutArray[i] - Min) / (Range);
		}
	}
}


void UHeightmapDeformersLibrary::ScalarOpArray(UPARAM(ref) TArray<float>& InOutArray, float Scale, EFloatAppendTypes AppendType /*= EFloatAppendTypes::Add*/)
{
	for (int32 i = 0; i < InOutArray.Num(); i++)
	{
		if (AppendType == EFloatAppendTypes::Add)
		{
			InOutArray[i] += Scale;
		}
		else if (AppendType == EFloatAppendTypes::Subtract)
		{
			InOutArray[i] -= Scale;
		}
		else if (AppendType == EFloatAppendTypes::Multiply)
		{
			InOutArray[i] *= Scale;
		}
		else if (AppendType == EFloatAppendTypes::Divide)
		{
			InOutArray[i] /= Scale;
		}
	}
}

TArray<float> UHeightmapDeformersLibrary::Conv_GreyScaleTexture2DToFloatArray(UTexture2D* InTexture)
{
	TArray<float> FloatArray;

	//sanity check for types we support atm
	if (InTexture->PlatformData->PixelFormat != PF_B8G8R8A8 &&
		InTexture->PlatformData->PixelFormat != PF_R8G8B8A8 && 
		InTexture->PlatformData->PixelFormat != PF_FloatRGBA)
	{
		UE_LOG(LogTemp, Warning, TEXT("Invalid float array conversion not yet supported for requested pixel format."));
		return FloatArray;
	}

	FloatArray.SetNum(InTexture->GetSizeX() * InTexture->GetSizeY());

	// Lock the texture so it can be read
	if (InTexture->PlatformData->PixelFormat == PF_FloatRGBA)
	{
		//Ensure settings for floaty rgba
		//InTexture->CompressionSettings = TextureCompressionSettings::TC_VectorDisplacementmap;
		InTexture->MipGenSettings = TextureMipGenSettings::TMGS_NoMipmaps;
		InTexture->SRGB = false;
		InTexture->UpdateResource();
		uint8* MipData = static_cast<uint8*>(InTexture->PlatformData->Mips[0].BulkData.Lock(LOCK_READ_ONLY));

		for (int i = 0; i < FloatArray.Num(); i++)
		{
			int MipPointer = i * 8;
			
			FFloat16Color GreyscaleValue;
			memcpy(&GreyscaleValue, &MipData[MipPointer], sizeof(FFloat16Color));

			//FloatArray[i] = (float)GreyscaleValue.r;

			FloatArray[i] = (float)GreyscaleValue.R;// / 65536.f;
			
			/// 65536.f;// / 255.f;	 //normalize it
		}
	}
	else
	{
		uint8* MipData = static_cast<uint8*>(InTexture->PlatformData->Mips[0].BulkData.Lock(LOCK_READ_ONLY));
		for (int32 i = 0; i < FloatArray.Num(); i++)
		{
			int32 MipPointer = i * 4;
			float GreyscaleValue = (MipData[MipPointer] + MipData[MipPointer + 1] + MipData[MipPointer + 2]) / 3.f;
			FloatArray[i] = GreyscaleValue / 255.f;	 //normalize it
		}
	}

	// Unlock the texture
	InTexture->PlatformData->Mips[0].BulkData.Unlock();

	return FloatArray;
}

UTexture2D* UHeightmapDeformersLibrary::Conv_GrayScaleFloatArrayToTexture2D(const TArray<float>& InFloatArray, const FVector2D InSize /*= FVector2D(0, 0)*/)
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

	UTexture2D* Pointer = SquareTextureSized(Size.X);
	CopyFloatArrayToTexture(InFloatArray, Pointer);

	return Pointer;
}


UTexture2D* UHeightmapDeformersLibrary::Conv_GrayScaleFloatArrayToHeightTexture2D(const TArray<float>& InFloatArray, const FVector2D InSize /*= FVector2D(0, 0)*/)
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

	UTexture2D* Pointer = SquareTextureSized(Size.X, EPixelFormat::PF_G16);
	CopyFloatArrayToTexture(InFloatArray, Pointer);

	return Pointer;
}

void UHeightmapDeformersLibrary::DeformTerrainByMask(TArray<float>& Patch, const TArray<float>& Mask, FTransform PatchTransform, FTransform MaskTransform, TFunction<float(float, float, float)> DeformAction, float Scale)
{
	//Assumption is square images only
	int32 PatchSize = (int32)FMath::Sqrt((float)Patch.Num());
	int32 MaskSize = (int32)FMath::Sqrt((float)Mask.Num());

	//Figure out scale by taking terrain scale vs mask scale. We assume they work in the same units
	//FVector DeltaScale = MaskTransform.GetScale3D() / TerrainTransform.GetScale3D();
	//FVector DeltaOrigin = (MaskTransform.GetLocation() * MaskTransform.GetScale3D()) - (PatchTransform.GetLocation()* PatchTransform.GetScale3D());

	FVector DeltaOrigin = PatchTransform.GetLocation() - MaskTransform.GetLocation();
	FVector DeltaScale = PatchTransform.GetScale3D() / MaskTransform.GetScale3D();

	for (float Y = 0.f; Y < PatchSize; Y++)
	{
		for (float X = 0.f; X < PatchSize; X++)
		{
			//Terrain index, read normally
			int32 PatchIndex = (PatchSize * Y) + X;

			if (PatchIndex >= Patch.Num())
			{
				UE_LOG(LogTemp, Log, TEXT("%d is out of bounds for Patch %d"), PatchIndex, Patch.Num())
				continue;
			}

			//MaskIndex, scale read, check validity
			float XSample = (X + DeltaOrigin.X) * DeltaScale.X;
			float YSample = FMath::RoundToInt((Y + DeltaOrigin.Y) * DeltaScale.Y) * MaskSize; //FMath::RoundToInt
			int32 MaskIndex = YSample + XSample;

			//This is expected to happen often, we skip these samples
			//Todo: allow unbounded/looped masking
			if (MaskIndex < 0 || Mask.Num() <= MaskIndex || XSample >= MaskSize || YSample < 0 || XSample < 0)
			{
				//UE_LOG(LogTemp, Log, TEXT("%d is out of bounds for Mask %d"), MaskIndex, Mask.Num())
				continue;
			}

			//Deform action returns actual desired value to give full control to deformer
			Patch[PatchIndex] = DeformAction(Patch[PatchIndex], Mask[MaskIndex], Scale);
		}
	}
}

void UHeightmapDeformersLibrary::DeformTerrainByMaskOp(TArray<float>& InOutTerrain, const TArray<float>& Mask, FTransform TerrainTransform, FTransform MaskTransform, EFloatAppendTypes MaskOp, float Scale /*= 1.f*/)
{
	//no op default
	TFunction<float(float, float, float)> MaskOpFunction = [](float TerrainPixel, float MaskPixel, float MaskScale) 
	{ 
		return TerrainPixel; 
	}; 
	
	if (MaskOp == EFloatAppendTypes::Add)
	{
		MaskOpFunction = [](float TerrainPixel, float MaskPixel, float MaskScale)
		{
			return TerrainPixel + (MaskPixel * MaskScale);
		};
	}
	else if (MaskOp == EFloatAppendTypes::Subtract)
	{
		MaskOpFunction = [](float TerrainPixel, float MaskPixel, float MaskScale)
		{
			return TerrainPixel - (MaskPixel * MaskScale);
		};
	}
	else if (MaskOp == EFloatAppendTypes::Multiply)
	{
		MaskOpFunction = [](float TerrainPixel, float MaskPixel, float MaskScale)
		{
			return TerrainPixel * (MaskPixel * MaskScale);
		};
	}
	else if (MaskOp == EFloatAppendTypes::Divide)
	{
		MaskOpFunction = [](float TerrainPixel, float MaskPixel, float MaskScale)
		{
			if (MaskPixel > 0.001f)	//avoid /0 problem
			{
				return TerrainPixel / (MaskPixel * MaskScale);
			}
			else
			{
				return TerrainPixel;
			}
		};
	}
	else if (MaskOp == EFloatAppendTypes::Override)
	{
		MaskOpFunction = [](float TerrainPixel, float MaskPixel, float MaskScale)
		{
			//If non-zero, override in mask area
			if (MaskPixel > 0.01f)
			{
				return MaskPixel * MaskScale;
			}
			else
			{
				return TerrainPixel;
			}
		};
	}

	DeformTerrainByMask(InOutTerrain, Mask, TerrainTransform, MaskTransform, MaskOpFunction, Scale);
}
