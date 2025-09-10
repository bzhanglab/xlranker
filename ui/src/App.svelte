<script lang="ts">
  import "./app.css";
  import Navbar from "$lib/components/Navbar.svelte";
  import { cn } from "$lib/utils"; // shadcn helper for conditional classes
  import NetworkUpload from "$lib/steps/NetworkUpload.svelte";
  import logo from "/src/assets/logo.svg";
  import Progress from "$lib/components/ui/progress/progress.svelte";

  let peptide_file: File[];
  let prog_value = 20;
  let currentStep = 1;
  let net_sep = "Tab";

  const steps = [
    { id: 1, label: "Upload Peptide Network" },
    { id: 2, label: "Omic Data Selection" },
    { id: 3, label: "Analysis" },
    { id: 4, label: "Results" },
  ];

  async function check_network() {
    const name = peptide_file.name;
    document.getElementById("status").textContent = "Starting task...";
    document.getElementById("result").textContent = "";

    const response = await fetch("/start_task", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name }),
    });

    const data = await response.json();
    currentTask = data.task_id;

    document.getElementById("status").textContent =
      "Task started. Checking progress...";
    pollTask();
  }
</script>

<Navbar />

<logo><img src={logo} alt="Logo" class="mx-auto mb-4 w-32" /></logo>

<div class="container mx-auto p-4 prose prose-lg max-w-3/4">
  <h1 class="text-center">XLRanker</h1>
  <div class="grid grid-cols-[200px_1fr] gap-6">
    <!-- Sidebar -->
    <aside class="flex flex-col space-y-1 text-sm">
      {#each steps as step}
        <div
          class={cn(
            "rounded-lg px-3 py-1 transition-colors",
            step.id === currentStep && "font-bold text-black",
            step.id < currentStep && "text-gray-700",
            step.id > currentStep && "text-gray-400"
          )}
        >
          {step.label}
        </div>
      {/each}
    </aside>
    <div>
      <div class="w-full grid space-x-4 text-sm px-4 py-4 rounded-lg bg-white">
        <NetworkUpload />
      </div>
      <!-- <Progress value={prog_value} max={100} id="progress" class="my-4" /> -->
    </div>
  </div>
</div>
