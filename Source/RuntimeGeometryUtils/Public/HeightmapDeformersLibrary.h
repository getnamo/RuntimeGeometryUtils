#pragma once

#include "Engine/Classes/Kismet/BlueprintFunctionLibrary.h"
#include "GeneratedMesh.h"
#include "HeightmapDeformersLibrary.generated.h"

USTRUCT()
struct FHeightAndGradient
{
	GENERATED_USTRUCT_BODY();

	UPROPERTY()
	float Height;

	UPROPERTY()
	float GradientX;

	UPROPERTY()
	float GradientY;
};

UENUM(BlueprintType)
enum class EFloatAppendTypes : uint8
{
	None,
	Add,
	Subtract,
	Multiply,
	Divide
};

USTRUCT(BlueprintType)
struct FHydroErosionParams
{
	GENERATED_USTRUCT_BODY();

	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	int32 Iterations;

	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	int32 Seed;

	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	int32 ErosionRadius;

	// At zero, water will instantly change direction to flow downhill. At 1, water will never change direction. 
	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	float Inertia;

	// Multiplier for how much sediment a droplet can carry
	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	float SedimentCapacityFactor;

	// Used to prevent carry capacity getting too close to zero on flatter terrain
	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	float MinSedimentCapacity;

	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	float ErodeSpeed;

	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	float DepositSpeed;

	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	float EvaporateSpeed;

	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	float Gravity;

	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	int32 MaxDropletLifetime;

	//Debug param: If true it will override output with 1.f and 0.f for deposit/erode
	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	bool bPreviewDropletErosionPaths;

	UPROPERTY(BlueprintReadWrite, EditAnywhere, Category = ErosionParams)
	bool bPreviewDropletDepositionPaths;

	FHydroErosionParams()
	{
		Iterations = 1;
		Seed = 1;
		ErosionRadius = 3;
		Inertia = 0.05f;
		SedimentCapacityFactor = 4;
		MinSedimentCapacity = 0.01f;
		ErodeSpeed = 0.3f;
		DepositSpeed = 0.3f;
		EvaporateSpeed = 0.01f;
		Gravity = 4;
		MaxDropletLifetime = 30;

		bPreviewDropletErosionPaths = false;
		bPreviewDropletDepositionPaths = false;
	}
};

/**
 * A BP Library of functions for applying deformations and erosions to heightmaps
 */
UCLASS(meta = (ScriptName = "GeneratedMeshDeformersLibrary"))
class RUNTIMEGEOMETRYUTILS_API UHeightmapDeformersLibrary : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()
public:

	//(Not really useful atm) PRNG via hashing. Seed will be modified each call
	UFUNCTION(BlueprintCallable) static
	int64 SplitMix64(int64& Seed);

	UFUNCTION(BlueprintCallable) static
	void PerlinDeformMap(UPARAM(ref) TArray<float>& InOutHeightmap, 
		float Magnitude = 1, 
		float Frequency = 1,
		FVector FrequencyShift = FVector(0, 0, 0), 
		int32 RandomSeed = 31337,
		int32 Octaves = 1,
		float OctaveFactor = 2.f, 
		bool bRidged = false);

	UFUNCTION(BlueprintCallable) static
	TArray<float> SquareFloatMapSized(int32 OneSideLength);

	UFUNCTION(BlueprintCallable) static
	UTexture2D* SquareTextureSized(int32 OneSideLength, EPixelFormat Format = PF_B8G8R8A8);

	UFUNCTION(BlueprintCallable) static
	void CopyFloatArrayToTexture(const TArray<float>& SrcData, UTexture2D* TargetTexture);

	/**
	* Erode a heightmap texture with hydraulic erosion via HydraulicErosionOnHeightMap. Convenience wrapper for texture conversion
	*/
	UFUNCTION(BlueprintCallable) static
	UTexture2D* HydraulicErosionOnHeightTexture(UTexture2D* InTexture, const FHydroErosionParams& Params);

	/** Erode a heightmap with particle simulation defined by a float array in place */
	UFUNCTION(BlueprintCallable) static
	void HydraulicErosionOnHeightMap(UPARAM(ref) TArray<float>& InOutHeightmap, const FHydroErosionParams& Params);

	UFUNCTION(BlueprintCallable) static
	UTexture2D* GenerateTransientCopy(UTexture2D* InTexture);

	// Not a super fan of the two below API-wise
	// consider unstable.

	// Used for e.g. append arrays
	UFUNCTION(BlueprintCallable) static
	void VectorOpArray(UPARAM(ref) TArray<float>& InOutArray, const TArray<float>& Other, EFloatAppendTypes AppendType = EFloatAppendTypes::Add, bool bNormalize = true);

	// Used for e.g. scaling arrays
	UFUNCTION(BlueprintCallable) static
	void ScalarOpArray(UPARAM(ref) TArray<float>& InOutArray, float Scale, EFloatAppendTypes AppendType = EFloatAppendTypes::Add);

	//We also need a sub-array append/scale with curve for blending

	//copies from https://github.com/getnamo/tensorflow-ue4/blob/master/Source/TensorFlow/Private/TensorFlowBlueprintLibrary.cpp
	UFUNCTION(BlueprintPure, meta = (DisplayName = "ToGrayScaleFloatArray (Texture2D)", BlueprintAutocast), Category = "Utilities|TensorFlow")
	static TArray<float> Conv_GreyScaleTexture2DToFloatArray(UTexture2D* InTexture);

	UFUNCTION(BlueprintPure, meta = (DisplayName = "ToTexture2D (Grayscale Array)", BlueprintAutocast), Category = "Utilities|TensorFlow")
	static UTexture2D* Conv_GrayScaleFloatArrayToTexture2D(const TArray<float>& InFloatArray, const FVector2D Size = FVector2D(0, 0));

	UFUNCTION(BlueprintPure, meta = (DisplayName = "ToHeightTexture2D (Grayscale Array)", BlueprintAutocast), Category = "Utilities|TensorFlow")
	static UTexture2D* Conv_GrayScaleFloatArrayToHeightTexture2D(const TArray<float>& InFloatArray, const FVector2D Size = FVector2D(0, 0));

	//UFUNCTION()


public:
	//C++ only
	/* Calculate gradient from map */
	static FHeightAndGradient CalculateHeightAndGradient(TArray<float>& Nodes, int32 MapSize, float PosX, float PosY);

	//Deform function signature: return value, terrain pix, mask pix
	static void DeformTerrainByMask(TArray<float>& InOutTerrain, const TArray<float>& Mask, FTransform MaskTransform, TFunction<float(float, float)> DeformAction);

	/**
	* Read from texture pointer and get height
	*/
	static float HeightAtPixel(float X, float Y, void* TexturePointer, int32 TextureHeight, int32 TextureWidth, int32 BytesPerPixel = 4);

	//Array needs to be square, Transform is relative to num squared unless scaled past 1.f
	static float HeightInSquareArray(const FTransform& SamplingTransform, TArray<float>& Array);

	static void HydraulicErosionOnHeightMapWithInterrupt(TArray<float>& InOutHeightmap, const FHydroErosionParams& Params, TFunction<bool()>InterruptHandler = nullptr);

protected:
	static void InitializeBrushIndices(int MapSize, int Radius, TArray<TArray<int32>>& ErosionBrushIndices, TArray<TArray<float>>& ErosionBrushWeights);
};