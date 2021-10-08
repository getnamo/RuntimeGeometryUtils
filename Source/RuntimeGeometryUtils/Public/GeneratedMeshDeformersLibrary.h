
#include "Engine/Classes/Kismet/BlueprintFunctionLibrary.h"
#include "GeneratedMesh.h"
#include "GeneratedMeshDeformersLibrary.generated.h"

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
 * A BP Library of functions for applying deformations to the vertices of a UGeneratedMesh
 */
UCLASS(meta = (ScriptName = "GeneratedMeshDeformersLibrary"))
class UGeneratedMeshDeformersLibrary : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()
public:

	/**
	 * Displace the mesh vertices using a 1D Sin wave.
	 * The displacement is based on a parameter t measured along the projection Axis.
	 * Once t is known, the new position is CurPos + Magnitude * Sin( (t*Frequency) + FrequencyShift ) * UpVector
	 */
	UFUNCTION(BlueprintCallable) static UPARAM(DisplayName = "Input Mesh") 
	UGeneratedMesh* DeformMeshAxisSinWave1D(UGeneratedMesh* Mesh, float Magnitude = 1, float Frequency = 1, float FrequencyShift = 0, FVector Axis = FVector(1,0,0), FVector UpVector = FVector(0,0,1));

	/**
	* Move vertices with specified radius/hardness to given location (e.g. g/move operator in blender)
	*/
	UFUNCTION(BlueprintCallable) static UPARAM(DisplayName = "Input Mesh")
	UGeneratedMesh* DeformMeshMove(UGeneratedMesh* Mesh, FVector FromWorldLocation = FVector(0, 0, 0), FVector ToWorldLocation = FVector(0, 0, 0), float Radius = 1, float Hardness = 1, float Magnitude = 1);


	/** 
	* Deform a mesh from a heightmap. Used to make procedural terrain from e.g. handcrafted heightmap
	*/
	UFUNCTION(BlueprintCallable) static UPARAM(DisplayName = "Input Mesh")
	UGeneratedMesh* DeformMeshHeightmap(UGeneratedMesh* Mesh, UTexture2D* Heightmap, float ZScale= 1, FVector DirectionNormal= FVector(0, 0, 1));


	/**
	 * Displace the mesh vertices using a 2D Sin wave.
	 */
	UFUNCTION(BlueprintCallable) static UPARAM(DisplayName = "Input Mesh") 
	UGeneratedMesh* DeformMeshAxisSinWaveRadial(UGeneratedMesh* Mesh, float Magnitude = 1, float Frequency = 1, float FrequencyShift = 0, FVector Axis = FVector(1,0,0));

	/**
	 * Displace the mesh vertices along their vertex normal directions using 3D Perlin Noise
	 */
	UFUNCTION(BlueprintCallable) static UPARAM(DisplayName = "Input Mesh")
	UGeneratedMesh* DeformMeshPerlinNoiseNormal(UGeneratedMesh* Mesh, float Magnitude = 1, float Frequency = 1, FVector FrequencyShift = FVector(0,0,0), int32 RandomSeed = 31337);

	/**
	 * Displace the mesh vertices along their vertex normal directions using 3D Perlin Noise, but only in narrow range of normal
	 */
	UFUNCTION(BlueprintCallable) static UPARAM(DisplayName = "Input Mesh")
	UGeneratedMesh* DeformMeshPerlinAlongAxisNormal(UGeneratedMesh* Mesh, float Magnitude = 1, float Frequency = 1, FVector FrequencyShift = FVector(0, 0, 0), int32 RandomSeed = 31337, FVector DirectionNormal = FVector(0, 0, 1));

	/**
	 * Displace the mesh vertices along their vertex normal directions using 3D Perlin Noise, but only above the given Z height (TODO: support any height cutoff vector)
	 */
	UFUNCTION(BlueprintCallable) static UPARAM(DisplayName = "Input Mesh")
	UGeneratedMesh* DeformMeshPerlinNormalPastHeight(UGeneratedMesh* Mesh, float Magnitude = 1, float Frequency = 1, FVector FrequencyShift = FVector(0, 0, 0), int32 RandomSeed = 31337, FVector HeightVector = FVector(0, 0, 1));

	/**
	 * Apply N iterations of explicit uniform Laplacian mesh smoothing to the vertex positions, with the given Alpha in range [0,1]. Clamps to max 100 iterations.
	 */
	UFUNCTION(BlueprintCallable) static UPARAM(DisplayName = "Input Mesh")
	UGeneratedMesh* SmoothMeshUniform(UGeneratedMesh* Mesh, float Alpha = 0.3, int32 Iterations = 1);

	//TODO: move these below to another bp library -> heightmap deform or 2.5 d deform ops

	UFUNCTION(BlueprintCallable) static
	void PerlinDeformMap(UPARAM(ref) TArray<float>& InOutHeightmap, float Magnitude = 1, float Frequency = 1, FVector FrequencyShift = FVector(0, 0, 0), int32 RandomSeed = 31337);

	UFUNCTION(BlueprintCallable) static
	TArray<float> SquareFloatMapSized(int32 OneSideLength);
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

	UFUNCTION(BlueprintCallable) static
	void AppendAndRescale(UPARAM(ref) TArray<float>& InOutArray, const TArray<float>& Other, bool bNormalize = true, EFloatAppendTypes AppendType = EFloatAppendTypes::Add);


	//copies from https://github.com/getnamo/tensorflow-ue4/blob/master/Source/TensorFlow/Private/TensorFlowBlueprintLibrary.cpp
	UFUNCTION(BlueprintPure, meta = (DisplayName = "ToGrayScaleFloatArray (Texture2D)", BlueprintAutocast), Category = "Utilities|TensorFlow")
	static TArray<float> Conv_GreyScaleTexture2DToFloatArray(UTexture2D* InTexture);

	UFUNCTION(BlueprintPure, meta = (DisplayName = "ToTexture2D (Grayscale Array)", BlueprintAutocast), Category = "Utilities|TensorFlow")
	static UTexture2D* Conv_GrayScaleFloatArrayToTexture2D(const TArray<float>& InFloatArray, const FVector2D Size = FVector2D(0, 0));

protected:
	/** 
	* Read from texture pointer and get height
	*/
	static float HeightAtPixel(float X, float Y, void* TexturePointer, int32 TextureHeight, int32 TextureWidth, int32 BytesPerPixel = 4);

	/* Calculate gradient from map */
	static FHeightAndGradient CalculateHeightAndGradient(TArray<float>& Nodes, int32 MapSize, float PosX, float PosY);
};