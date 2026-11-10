#!/bin/bash

fixed="/mnt/c/Users/78457/Documents/MNI152/MNI152_T1_2mm.nii.gz"

input_dir="/mnt/c/Users/78457/Documents/MNI152/IXI-T1"
output_dir="/mnt/c/Users/78457/Documents/MNI152/IXI-T1/registration_batch"

mkdir -p "${output_dir}/transforms"
mkdir -p "${output_dir}/mni_images"

for moving in "${input_dir}"/*T1.nii.gz
do
    filename=$(basename "${moving}")
    subject="${filename%%-*}"

    echo "========================================"
    echo "Processing ${subject}"
    echo "Input: ${filename}"
    echo "========================================"

    output_prefix="${output_dir}/transforms/${subject}_to_MNI_"
    output_image="${output_dir}/mni_images/${subject}_MNI_T1_2mm.nii.gz"

    if [ -f "${output_image}" ]; then
        echo "${subject} already finished. Skipping."
        continue
    fi

    echo "Registering ${subject}..."

    antsRegistration \
      --dimensionality 3 \
      --float 0 \
      --output "${output_prefix}" \
      --interpolation Linear \
      --use-histogram-matching 0 \
      --winsorize-image-intensities "[0.005,0.995]" \
      --initial-moving-transform "[${fixed},${moving},1]" \
      --transform Rigid[0.1] \
      --metric "MI[${fixed},${moving},1,32,Regular,0.25]" \
      --convergence "[1000x500x250x100,1e-6,10]" \
      --shrink-factors 8x4x2x1 \
      --smoothing-sigmas 3x2x1x0vox \
      --transform Affine[0.1] \
      --metric "MI[${fixed},${moving},1,32,Regular,0.25]" \
      --convergence "[1000x500x250x100,1e-6,10]" \
      --shrink-factors 8x4x2x1 \
      --smoothing-sigmas 3x2x1x0vox \
      --transform SyN[0.1,3,0] \
      --metric "CC[${fixed},${moving},1,4]" \
      --convergence "[100x70x50x20,1e-6,10]" \
      --shrink-factors 8x4x2x1 \
      --smoothing-sigmas 3x2x1x0vox

    if [ $? -ne 0 ]; then
        echo "ERROR: Registration failed for ${subject}."
        continue
    fi

    echo "Applying transforms..."

    antsApplyTransforms \
      -d 3 \
      -i "${moving}" \
      -r "${fixed}" \
      -o "${output_image}" \
      -n Linear \
      -t "${output_prefix}1Warp.nii.gz" \
      -t "${output_prefix}0GenericAffine.mat"

    if [ $? -ne 0 ]; then
        echo "ERROR: Transform failed for ${subject}."
        continue
    fi

    echo "${subject} finished successfully."
done

echo "========================================"
echo "All subjects processed."
echo "========================================"
