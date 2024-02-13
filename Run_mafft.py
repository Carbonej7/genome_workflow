import sys
import subprocess


# This will Run MAFFT alignment and remove Reference sequence.

def run_mafft(input_file, ref_sequence, output_file):
    cmd = ["mafft", "--6merpair", "--thread", "4", "--keeplength"]
    if ref_sequence:
        cmd.extend(["--addfragments", ref_sequence])
    cmd.append(input_file)

    try:
        with open(output_file, "w") as f:
            subprocess.run(cmd, check=True, stdout=f)
        print("MAFFT alignment completed successfully.")

        if ref_sequence:
            with open(output_file, "r") as f:
                lines = f.readlines()
            with open(output_file, "w") as f:
                f.writelines(lines[1:])
    except subprocess.CalledProcessError as e:
        print(f"Error: MAFFT alignment failed with exit code {e.returncode}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python Run_mafft.py input_sequences input_ref output_sequences")
    input_file = sys.argv[1]
    ref_sequence = sys.argv[2]
    output_file = sys.argv[3]

    run_mafft(ref_sequence, input_file, output_file)
