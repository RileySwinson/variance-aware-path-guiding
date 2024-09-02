import subprocess

result = subprocess.run(["mtsutil", "dscompare"])

print(result.stdout)