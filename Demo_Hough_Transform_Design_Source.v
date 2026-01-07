`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 01/07/2026 09:47:43 AM
// Design Name: 
// Module Name: Demo_Hough_Transform_Design_Source
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module Demo_Hough_Transform_Design_Source #(

    parameter IMG_WIDTH         = 640,
    parameter IMG_HEIGHT        = 480,
    parameter IO_DATA_WIDTH     = 24,
    parameter ROI_Y_START       = 200,
    parameter ROI_Y_END         = 479,
    parameter THETA_STEPS       = 180,
    parameter RHO_MAX           = 800,
    parameter ACCUM_WIDTH       = 8, 
    parameter THETA_DIVIDER     = 90
    )
    (
    
    input  wire                       aresetn,
    input  wire                       ap_clk,
    input  wire [IO_DATA_WIDTH-1:0]   s_axis_tdata,
    input  wire                       s_axis_tvalid,
    output wire                       s_axis_tready,
    input  wire                       s_axis_tlast,
    input  wire                       s_axis_tuser,
    output wire [IO_DATA_WIDTH-1:0]   m_axis_tdata,
    output wire                       m_axis_tvalid,
    input  wire                       m_axis_tready,
    output wire                       m_axis_tlast,
    output wire                       m_axis_tuser    
    
    );
    
localparam FP_BITS             = 8;
localparam CALC_WIDTH          = 24;
localparam LUT_WIDTH           = 16;
localparam FRAME_SIZE          = IMG_WIDTH * IMG_HEIGHT;
localparam FRAME_ADDR_WIDTH    = $clog2(FRAME_SIZE);
localparam FIFO_X_WIDTH        = $clog2(IMG_WIDTH)+1;
localparam FIFO_Y_WIDTH        = $clog2(IMG_HEIGHT)+1;
localparam FIFO_COORD_WIDTH    = FIFO_X_WIDTH + FIFO_Y_WIDTH;
localparam FIFO_DEPTH          = 65536/4;
localparam LINE_THICKNESS      = 3;
localparam Y_SKIP_INTERVAL     = 5;
localparam COS_SIN_FRAC_BITS   = 11;
localparam RHO_DIAMETER        = 2 * RHO_MAX;
localparam ACCUM_SIZE          = THETA_STEPS * RHO_DIAMETER;
localparam ACCUM_ADDR_WIDTH    = $clog2(ACCUM_SIZE);
localparam signed [31:0]       CENTER_Y_RHO_TRUE = $signed(-((IMG_HEIGHT-8) / 2));
localparam [$clog2(RHO_DIAMETER)-1:0] CENTER_Y_RHO_IDX = CENTER_Y_RHO_TRUE + RHO_MAX;
localparam signed [64:0]       ROUNDING_CONSTANT = (1 << (COS_SIN_FRAC_BITS - 1));
localparam SHIFT_CONVERSION    = COS_SIN_FRAC_BITS - FP_BITS;

localparam S_IDLE              = 5'd0;
localparam S_RESET_ACCUM       = 5'd1;
localparam S_RECEIVE_FRAME     = 5'd2;
localparam S_HOUGH_ACCUM       = 5'd3;
localparam S_FIND_PEAKS        = 5'd4;
localparam S_CALCULATE_LINES   = 5'd5;
localparam S_DRAW_FRAME        = 5'd6;

localparam PIPE_FETCH_COORD    = 5'd0; 
localparam PIPE_LATCH_XY_AND_INIT_THETA = 5'd1; 
localparam PIPE_CALC_MUL_L0    = 5'd2; 
localparam PIPE_CALC_SUM_SHIFT = 5'd3; 
localparam PIPE_CALC_ADDR      = 5'd4; 
localparam PIPE_READ_HOUGH_ACCUM = 5'd5; 
localparam PIPE_MODIFY_ACCUM   = 5'd6;
localparam PIPE_WRITE_HOUGH_ACCUM = 5'd7; 
localparam PIPE_LOOP_THETA_ITER = 5'd8; 

localparam NMS_THETA_WINDOW    = 20;
localparam NMS_RHO_WINDOW      = 20;
localparam THETA_LANE_MIN      = 1;
localparam THETA_LANE_MAX      = 179;
localparam THETA_STAB_THRESHOLD = 30;
localparam RHO_STAB_THRESHOLD   = 150;
localparam BLACK               = 24'h000000;
localparam WHITE               = 24'hFFFFFF;
localparam LANE_COLOR          = 24'h0000ff;
localparam MIN_PEAK_THRESHOLD  = 25; 
localparam TOTAL_DRAW_SHIFT    = 2 * COS_SIN_FRAC_BITS;




reg [4:0] state, next_state;
reg       frame_done_flag;
reg       reset_done_flag;
reg       accum_done_flag;
reg       peaks_done_flag;
reg       calc_done_flag;
reg       draw_done_flag;

reg [FIFO_X_WIDTH-1:0]         x_in_cnt;
reg [FIFO_Y_WIDTH-1:0]         y_in_cnt;
reg [FIFO_X_WIDTH-1:0]         x_out_cnt;
reg [FIFO_Y_WIDTH-1:0]         y_out_cnt;

reg [3:0]                      accum_pipe_state;
reg                            accum_processing;
(* mark_debug = "true" *) reg [FIFO_COORD_WIDTH-1:0]     fifo_coord_reg;

reg signed [FIFO_X_WIDTH-1:0]  hc_x_coord_aligned;
reg signed [FIFO_Y_WIDTH-1:0]  hc_y_coord_aligned;

reg [$clog2(THETA_STEPS)-1:0]  accum_theta_idx;          
reg [$clog2(THETA_STEPS)-1:0]  accum_theta_idx_aligned;  

reg [ACCUM_ADDR_WIDTH-1:0]     accum_rd_addr_reg;
reg [ACCUM_ADDR_WIDTH-1:0]     accum_wr_addr_reg;
reg [ACCUM_WIDTH-1:0]          accum_read_data_reg;

reg signed [47:0]              hc_x_cos_prod_reg;
reg signed [47:0]              hc_y_sin_prod_reg;
reg signed [64:0]              hc_rho_sum_shifted_reg;
reg [$clog2(THETA_STEPS)-1:0]  accum_theta_idx_p1;
reg [$clog2(THETA_STEPS)-1:0]  accum_theta_idx_p2;

reg [$clog2(THETA_STEPS)-1:0]  find_theta_cnt;
reg [$clog2(RHO_DIAMETER)-1:0] find_rho_cnt;
reg [ACCUM_WIDTH-1:0]          hough_accum_dout_b_reg;
reg [$clog2(THETA_STEPS)-1:0]  find_theta_cnt_reg;
reg [$clog2(RHO_DIAMETER)-1:0] find_rho_cnt_reg;
reg [ACCUM_ADDR_WIDTH-1:0]     find_peaks_addr_reg;
reg [ACCUM_ADDR_WIDTH-1:0]     reset_addr_reg;

reg [ACCUM_WIDTH-1:0]          peak1_val_nms;
reg [ACCUM_WIDTH-1:0]          peak2_val_nms; 
reg [$clog2(THETA_STEPS)-1:0]  peak1_theta_nms;
reg [$clog2(RHO_DIAMETER)-1:0] peak1_rho_idx_nms;
reg [$clog2(THETA_STEPS)-1:0]  peak2_theta_nms;
reg [$clog2(RHO_DIAMETER)-1:0] peak2_rho_idx_nms;

reg [3:0]                      calc_step_cnt;
reg [$clog2(THETA_STEPS)-1:0]  sorted_theta_left;
reg [$clog2(RHO_DIAMETER)-1:0] sorted_rho_idx_left;
reg [$clog2(THETA_STEPS)-1:0]  sorted_theta_right;
reg [$clog2(RHO_DIAMETER)-1:0] sorted_rho_idx_right;

reg [$clog2(THETA_STEPS)-1:0]  current_cand_theta_left;
reg [$clog2(RHO_DIAMETER)-1:0] current_cand_rho_idx_left;
reg [$clog2(THETA_STEPS)-1:0]  current_cand_theta_right;
reg [$clog2(RHO_DIAMETER)-1:0] current_cand_rho_idx_right;

reg [$clog2(THETA_STEPS)-1:0]  final_theta_left;
reg [$clog2(RHO_DIAMETER)-1:0] final_rho_idx_left;
reg [$clog2(THETA_STEPS)-1:0]  final_theta_right;
reg [$clog2(RHO_DIAMETER)-1:0] final_rho_idx_right;

reg [$clog2(THETA_STEPS)-1:0]  prev_final_theta_left_r;
reg [$clog2(RHO_DIAMETER)-1:0] prev_final_rho_idx_left_r;
reg                             prev_final_left_valid_r;
reg [$clog2(THETA_STEPS)-1:0]  prev_final_theta_right_r;
reg [$clog2(RHO_DIAMETER)-1:0] prev_final_rho_idx_right_r;
reg                             prev_final_right_valid_r;

reg signed [47:0]              wide_rho_true_left;
reg signed [47:0]              wide_rho_true_right;

// Output wires from BRAMs
 wire signed [15:0] cos_val_bram_out_accum_pipe;
 wire signed [15:0] sin_val_bram_out_accum_pipe;

 wire signed [15:0] sec_val_bram_out_calc_left;
 wire signed [15:0] neg_tan_val_bram_out_calc_left;
 wire signed [15:0] sec_val_bram_out_calc_right;
 wire signed [15:0] neg_tan_val_bram_out_calc_right;





wire signed [15:0] sin_val_bram_out_draw_left;
wire signed [15:0] sin_val_bram_out_draw_right;

reg signed [47:0] wide_sec_lut_left_aligned;
reg signed [47:0] wide_neg_tan_lut_left_aligned;
reg signed [47:0] wide_sec_lut_right_aligned;
reg signed [47:0] wide_neg_tan_lut_right_aligned;
reg signed [47:0] sin_val_draw_left_aligned;
reg signed [47:0] sin_val_draw_right_aligned;


reg signed [FIFO_X_WIDTH:0]    x_lane_left;
reg signed [FIFO_X_WIDTH:0]    x_lane_right;

reg       write_buffer_idx;
reg       read_buffer_idx;
reg       read_buffer_valid;

reg signed [47:0] draw_y_sin_prod_left;
reg signed [47:0] draw_y_sin_prod_right;
reg signed [64:0] draw_numerator_left;
reg signed [64:0] draw_numerator_right;
reg signed [64:0] draw_product_left;
reg signed [64:0] draw_product_right;
reg signed [47:0] left_term1;
reg signed [47:0] left_term2;
reg signed [47:0] right_term1;
reg signed [47:0] right_term2;
reg signed [48:0] left_product;
reg signed [48:0] right_product;

wire fb_we_a;
wire hough_accum_we_a;
wire binary_pixel_in;

wire binary_pixel_out_b;
wire fifo_wr_en;
(* mark_debug = "true" *) reg fifo_wr_en_db;
wire fifo_rd_en;
(* mark_debug = "true" *) reg fifo_rd_en_db;
wire fifo_full;
(* mark_debug = "true" *) reg fifo_full_db;
wire fifo_empty;
(* mark_debug = "true" *) reg fifo_empty_db;


wire [FRAME_ADDR_WIDTH:0]      fb_addr_a;
wire [FRAME_ADDR_WIDTH:0]      fb_addr_b;
wire [ACCUM_ADDR_WIDTH-1:0]    hough_accum_addr_a;
wire [ACCUM_ADDR_WIDTH-1:0]    hough_accum_addr_b;
wire [ACCUM_WIDTH-1:0]         hough_accum_din_a;
wire [ACCUM_WIDTH-1:0]         hough_accum_dout_b;
wire [FIFO_COORD_WIDTH-1:0]    fifo_rd_data;


always @(*)begin

fifo_wr_en_db = fifo_wr_en;
fifo_full_db = fifo_full;

fifo_empty_db = fifo_empty;
fifo_rd_en_db = fifo_rd_en;
end

cos_mem_gen_0 i_cos_lut_accum (
    .clka(ap_clk),
    .addra(accum_theta_idx),
    .douta(cos_val_bram_out_accum_pipe)
);
sin_mem_gen_0 i_sin_lut_accum (
    .clka(ap_clk),
    .addra(accum_theta_idx),
    .douta(sin_val_bram_out_accum_pipe)
);



sec_mem_gen_0 i_sec_lut_calc_left (
    .clka(ap_clk),
    .addra(final_theta_left),
    .douta(sec_val_bram_out_calc_left)
);
neg_tan_mem_gen_0 i_neg_tan_lut_calc_left (
    .clka(ap_clk),
    .addra(final_theta_left),
    .douta(neg_tan_val_bram_out_calc_left)
);

sec_mem_gen_0 i_sec_lut_calc_right (
    .clka(ap_clk),
    .addra(final_theta_right),
    .douta(sec_val_bram_out_calc_right)
);
neg_tan_mem_gen_0 i_neg_tan_lut_calc_right (
    .clka(ap_clk),
    .addra(final_theta_right),
    .douta(neg_tan_val_bram_out_calc_right)
);


sin_mem_gen_0 i_sin_lut_draw_left (
    .clka(ap_clk),
    .addra(final_theta_left),
    .douta(sin_val_bram_out_draw_left)
);
sin_mem_gen_0 i_sin_lut_draw_right (
    .clka(ap_clk),
    .addra(final_theta_right),
    .douta(sin_val_bram_out_draw_right)
);



frame_buffer_bram frame_buffer_inst (
    .clka(ap_clk),
    .wea(fb_we_a),
    .addra(fb_addr_a),
    .dina(binary_pixel_in),
    .clkb(ap_clk),
    .addrb(fb_addr_b),
    .doutb(binary_pixel_out_b)
);

hough_accum_bram hough_accum_inst (
    .clka(ap_clk),
    .wea(hough_accum_we_a),
    .addra(hough_accum_addr_a),
    .dina(hough_accum_din_a),
    .clkb(ap_clk),
    .addrb(hough_accum_addr_b),
    .doutb(hough_accum_dout_b)
);

Demo_fifo_Design_Source #(.DATA_WIDTH(FIFO_COORD_WIDTH), .DEPTH(FIFO_DEPTH)) coord_fifo (
    .clk(ap_clk),
    .resetn(aresetn),
    .wr_data({y_in_cnt, x_in_cnt}),
    .wr_en(fifo_wr_en),
    .full(fifo_full),
    .rd_data(fifo_rd_data),
    .rd_en(fifo_rd_en),
    .empty(fifo_empty)
);

assign binary_pixel_in = s_axis_tdata[0];
assign fb_addr_a       = {write_buffer_idx, (y_in_cnt * IMG_WIDTH + x_in_cnt)};
assign fb_addr_b       = {read_buffer_idx,  (y_out_cnt * IMG_WIDTH + x_out_cnt)};


wire signed [FIFO_X_WIDTH-1:0] hc_x_coord_from_fifo = fifo_coord_reg[FIFO_X_WIDTH-1:0];
wire signed [FIFO_Y_WIDTH-1:0] hc_y_coord_from_fifo = fifo_coord_reg[FIFO_COORD_WIDTH-1:FIFO_X_WIDTH];


wire signed [47:0]             x_op_aligned = $signed(hc_x_coord_aligned);
wire signed [47:0]             y_op_aligned = $signed(hc_y_coord_aligned);


wire signed [47:0]             cos_op_for_accum = $signed(cos_val_bram_out_accum_pipe);
wire signed [47:0]             sin_op_for_accum = $signed(sin_val_bram_out_accum_pipe);


wire signed [47:0]             hc_x_cos_prod         = x_op_aligned * cos_op_for_accum;
wire signed [47:0]             hc_y_sin_prod         = y_op_aligned * sin_op_for_accum;

wire signed [64:0]             hc_rho_sum_unrounded  = {{16{hc_x_cos_prod_reg[47]}}, hc_x_cos_prod_reg} +
                                                        {{16{hc_y_sin_prod_reg[47]}}, hc_y_sin_prod_reg};
wire signed [64:0]             hc_rho_sum_rounded    = hc_rho_sum_unrounded + ROUNDING_CONSTANT;
wire signed [64:0]             hc_rho_sum_shifted    = hc_rho_sum_rounded >>> COS_SIN_FRAC_BITS;
wire signed [31:0]             hc_rho_val            = hc_rho_sum_shifted_reg;
wire signed [31:0]             hc_temp_rho_idx_signed = hc_rho_val + RHO_MAX;
wire [$clog2(RHO_DIAMETER)-1:0] hc_temp_rho_idx = (hc_temp_rho_idx_signed < 0) ? 0 :
                                                   (hc_temp_rho_idx_signed >= RHO_DIAMETER) ? (RHO_DIAMETER-1) :
                                                   hc_temp_rho_idx_signed;
wire [ACCUM_ADDR_WIDTH-1:0]    calculated_addr       = accum_theta_idx_p2 * RHO_DIAMETER + hc_temp_rho_idx;

wire signed [$clog2(THETA_STEPS):0]  theta_diff_stab_l = $signed(current_cand_theta_left) - $signed(prev_final_theta_left_r);
wire signed [$clog2(RHO_DIAMETER):0] rho_diff_stab_l = $signed(current_cand_rho_idx_left) - $signed(prev_final_rho_idx_left_r);
wire [$clog2(THETA_STEPS)-1:0]       abs_theta_diff_stab_l = (theta_diff_stab_l < 0) ? -theta_diff_stab_l : theta_diff_stab_l;

wire [$clog2(RHO_DIAMETER)-1:0]      abs_rho_diff_stab_l = (rho_diff_stab_l < 0) ? -rho_diff_stab_l : rho_diff_stab_l;
wire is_left_lane_stable = (prev_final_left_valid_r == 1) &&
                           ((abs_theta_diff_stab_l <= THETA_STAB_THRESHOLD) &&
                            (abs_rho_diff_stab_l <= RHO_STAB_THRESHOLD));

wire signed [$clog2(THETA_STEPS):0]  theta_diff_stab_r = $signed(current_cand_theta_right) - $signed(prev_final_theta_right_r);
wire signed [$clog2(RHO_DIAMETER):0] rho_diff_stab_r = $signed(current_cand_rho_idx_right) - $signed(prev_final_rho_idx_right_r);
wire [$clog2(THETA_STEPS)-1:0]       abs_theta_diff_stab_r = (theta_diff_stab_r < 0) ? -theta_diff_stab_r : theta_diff_stab_r;

wire [$clog2(RHO_DIAMETER)-1:0]      abs_rho_diff_stab_r = (rho_diff_stab_r < 0) ? -rho_diff_stab_r : rho_diff_stab_r;
wire is_right_lane_stable = (prev_final_right_valid_r == 1) &&
                            ((abs_theta_diff_stab_r <= THETA_STAB_THRESHOLD) &&
                             (abs_rho_diff_stab_r <= RHO_STAB_THRESHOLD));


always @(*) begin
    next_state = state;
    case(state)
        S_IDLE:              next_state = S_RESET_ACCUM;
        S_RESET_ACCUM:       if (reset_done_flag)    next_state = S_RECEIVE_FRAME;
        S_RECEIVE_FRAME:     if (frame_done_flag)    next_state = S_HOUGH_ACCUM;
        S_HOUGH_ACCUM:       if (accum_done_flag)    next_state = S_FIND_PEAKS;
        S_FIND_PEAKS:        if (peaks_done_flag)    next_state = S_CALCULATE_LINES;
        S_CALCULATE_LINES:   if (calc_done_flag)     next_state = S_DRAW_FRAME;
        S_DRAW_FRAME:        if (draw_done_flag)     next_state = S_IDLE;
        default:             next_state = S_IDLE;
    endcase
end

assign s_axis_tready = (state == S_RECEIVE_FRAME) && !fifo_full;
assign m_axis_tvalid = (state == S_DRAW_FRAME) && read_buffer_valid;
assign m_axis_tuser  = (x_out_cnt == 0) && (y_out_cnt == 0);
assign m_axis_tlast  = (x_out_cnt == IMG_WIDTH-1);

assign fb_we_a       = (state == S_RECEIVE_FRAME) && (s_axis_tvalid && s_axis_tready);
assign fifo_wr_en    = (state == S_RECEIVE_FRAME) && (s_axis_tvalid && s_axis_tready) &&
                       (y_in_cnt >= ROI_Y_START && y_in_cnt <= ROI_Y_END) &&
                       ((y_in_cnt % Y_SKIP_INTERVAL) == 0) && binary_pixel_in && !fifo_full;
assign fifo_rd_en    = (state == S_HOUGH_ACCUM) && (accum_pipe_state == PIPE_FETCH_COORD) &&
                       !accum_processing && !fifo_empty;

assign hough_accum_we_a = (state == S_RESET_ACCUM) ? 1'b1 :
                          (state == S_HOUGH_ACCUM && accum_pipe_state == PIPE_WRITE_HOUGH_ACCUM) ? 1'b1 : 1'b0;
assign hough_accum_addr_a = (state == S_RESET_ACCUM) ? reset_addr_reg : accum_wr_addr_reg;
assign hough_accum_addr_b = (state == S_FIND_PEAKS) ? find_peaks_addr_reg : accum_rd_addr_reg;
assign hough_accum_din_a = (state == S_RESET_ACCUM) ? {ACCUM_WIDTH{1'b0}} : (hough_accum_dout_b + 1);

wire is_on_right_lane = prev_final_right_valid_r &&
                        (x_lane_right >= 0) && (x_lane_right < IMG_WIDTH) &&
                        (y_out_cnt >= ROI_Y_START) &&
                        ($signed(x_out_cnt) >= (x_lane_right - LINE_THICKNESS)) &&
                        ($signed(x_out_cnt) <= (x_lane_right + LINE_THICKNESS));

wire is_on_left_lane  = prev_final_left_valid_r &&
                        (x_lane_left >= 0) && (x_lane_left < IMG_WIDTH) &&
                        (y_out_cnt >= ROI_Y_START) &&
                        ($signed(x_out_cnt) >= (x_lane_left - LINE_THICKNESS)) &&
                        ($signed(x_out_cnt) <= (x_lane_left + LINE_THICKNESS));

assign m_axis_tdata = (is_on_right_lane || is_on_left_lane) ? LANE_COLOR :
                      (binary_pixel_out_b) ? WHITE : BLACK;


always @(posedge ap_clk or negedge aresetn) begin
    if (!aresetn) begin
        state                      <= S_IDLE;
        x_in_cnt                   <= 0;
        y_in_cnt                   <= 0;
        x_out_cnt                  <= 0;
        y_out_cnt                  <= 0;
        frame_done_flag            <= 0;
        reset_done_flag            <= 0;
        accum_done_flag            <= 0;
        peaks_done_flag            <= 0;
        calc_done_flag             <= 0;
        draw_done_flag             <= 0;
        write_buffer_idx           <= 0;
        read_buffer_idx            <= 1;
        read_buffer_valid          <= 0;
        fifo_coord_reg             <= 0;
        hc_x_coord_aligned         <= 0;
        hc_y_coord_aligned         <= 0;
        accum_theta_idx            <= 0;
        accum_theta_idx_aligned    <= 0;
        accum_pipe_state           <= PIPE_FETCH_COORD;
        accum_read_data_reg        <= 0;
        accum_processing           <= 1'b0;
        find_theta_cnt             <= 0;
        find_rho_cnt               <= 0;
        hough_accum_dout_b_reg     <= 0;
        find_theta_cnt_reg         <= 0;
        find_rho_cnt_reg           <= 0;
        peak1_val_nms              <= 0;
        peak2_val_nms              <= 0;
        peak1_theta_nms            <= 0;
        peak2_theta_nms            <= 0;
        peak1_rho_idx_nms          <= 0;
        peak2_rho_idx_nms          <= 0;
        hc_x_cos_prod_reg          <= 0;
        hc_y_sin_prod_reg          <= 0;
        hc_rho_sum_shifted_reg     <= 0;
        accum_theta_idx_p1         <= 0;
        accum_theta_idx_p2         <= 0;
        accum_rd_addr_reg          <= 0;
        accum_wr_addr_reg          <= 0;
        reset_addr_reg             <= 0;
        find_peaks_addr_reg        <= 0;
        prev_final_theta_left_r    <= 0;
        prev_final_rho_idx_left_r  <= 0;
        prev_final_left_valid_r    <= 0;
        prev_final_theta_right_r   <= 0;
        prev_final_rho_idx_right_r <= 0;
        prev_final_right_valid_r   <= 0;
        sorted_theta_left          <= 0;
        sorted_rho_idx_left        <= 0;
        sorted_theta_right         <= 0;
        sorted_rho_idx_right       <= 0;
        current_cand_theta_left    <= 0;
        current_cand_theta_right   <= 0;
        current_cand_rho_idx_left  <= 0;
        current_cand_rho_idx_right <= 0;
        final_theta_left           <= 0;
        final_theta_right          <= 0;
        final_rho_idx_left         <= 0;
        final_rho_idx_right        <= 0;
        calc_step_cnt              <= 0;
        wide_rho_true_left         <= 0;
        wide_rho_true_right        <= 0;
        wide_sec_lut_left_aligned      <= 0;
        wide_neg_tan_lut_left_aligned  <= 0;
        wide_sec_lut_right_aligned     <= 0;
        wide_neg_tan_lut_right_aligned <= 0;
        sin_val_draw_left_aligned      <= 0;
        sin_val_draw_right_aligned     <= 0;
        x_lane_left                <= 0;
        x_lane_right               <= 0;
        draw_y_sin_prod_left       <= 0;
        draw_y_sin_prod_right      <= 0;
        draw_numerator_left        <= 0;
        draw_numerator_right       <= 0;
        draw_product_left          <= 0;
        draw_product_right         <= 0;
        left_term1                 <= 0;
        left_term2                 <= 0;
        right_term1                <= 0;
        right_term2                <= 0;
        left_product               <= 0;
        right_product              <= 0;
    end else begin
        frame_done_flag     <= 0;
        reset_done_flag     <= 0;
        accum_done_flag     <= 0;
        peaks_done_flag     <= 0;
        calc_done_flag      <= 0;
        draw_done_flag      <= 0;
        state               <= next_state;

        hough_accum_dout_b_reg <= hough_accum_dout_b;
        find_theta_cnt_reg     <= find_theta_cnt;
        find_rho_cnt_reg       <= find_rho_cnt;

        if (calc_done_flag) begin
            prev_final_theta_left_r    <= final_theta_left;
            prev_final_rho_idx_left_r  <= final_rho_idx_left;
            prev_final_left_valid_r    <= (peak2_val_nms >= MIN_PEAK_THRESHOLD);

            prev_final_theta_right_r   <= final_theta_right;
            prev_final_rho_idx_right_r <= final_rho_idx_right;
            prev_final_right_valid_r   <= (peak1_val_nms >= MIN_PEAK_THRESHOLD);
        end

        case(state)
            S_IDLE: begin
                find_theta_cnt <= 0;
                find_rho_cnt   <= 0;
            end

            S_RESET_ACCUM: begin
                reset_addr_reg <= find_theta_cnt * RHO_DIAMETER + find_rho_cnt;
                if (find_rho_cnt == RHO_DIAMETER - 1) begin
                    find_rho_cnt <= 0;
                    if (find_theta_cnt == THETA_STEPS - 1) begin
                        reset_done_flag <= 1'b1;
                        x_in_cnt        <= 0;
                        y_in_cnt        <= 0;
                    end else begin
                        find_theta_cnt <= find_theta_cnt + 1;
                    end
                end else begin
                    find_rho_cnt <= find_rho_cnt + 1;
                end
            end

            S_RECEIVE_FRAME: begin
                if (s_axis_tvalid && s_axis_tready) begin
                    if (x_in_cnt == IMG_WIDTH - 1 && y_in_cnt == IMG_HEIGHT - 1) begin
                        frame_done_flag   <= 1'b1;
                        read_buffer_idx   <= write_buffer_idx;
                        write_buffer_idx  <= ~write_buffer_idx;
                        read_buffer_valid <= 1'b1;
                        x_in_cnt          <= 0;
                        y_in_cnt          <= 0;
                    end else begin
                        if (x_in_cnt == IMG_WIDTH - 1) begin
                            x_in_cnt <= 0;
                            y_in_cnt <= y_in_cnt + 1;
                        end else begin
                            x_in_cnt <= x_in_cnt + 1;
                        end
                    end
                end
            end

            S_HOUGH_ACCUM: begin
                case (accum_pipe_state)
                    PIPE_FETCH_COORD: begin
                        if (!accum_processing && !fifo_empty) begin
                            accum_pipe_state <= PIPE_LATCH_XY_AND_INIT_THETA;
                        end else if (!accum_processing && fifo_empty) begin
                            accum_done_flag      <= 1'b1;
                            find_theta_cnt       <= 0;
                            find_rho_cnt         <= 0;
                            peak1_val_nms        <= 0;
                            peak2_val_nms        <= 0;
                            peak1_theta_nms      <= 0;
                            peak2_theta_nms      <= 0;
                            peak1_rho_idx_nms    <= 0;
                            peak2_rho_idx_nms    <= 0;
                        end
                    end
                    PIPE_LATCH_XY_AND_INIT_THETA: begin 
                        fifo_coord_reg   <= fifo_rd_data;
                        hc_x_coord_aligned <= hc_x_coord_from_fifo;
                        hc_y_coord_aligned <= hc_y_coord_from_fifo;
                        
                        accum_processing <= 1'b1;
                        accum_theta_idx  <= 0; 
                        accum_theta_idx_aligned <= 0; 
                        accum_pipe_state <= PIPE_CALC_MUL_L0; 
                    end
                    PIPE_CALC_MUL_L0: begin 
                        hc_x_cos_prod_reg  <= x_op_aligned * cos_op_for_accum; 
                        hc_y_sin_prod_reg  <= y_op_aligned * sin_op_for_accum; 
                        accum_theta_idx_p1 <= accum_theta_idx_aligned;
                        

                        if (accum_theta_idx_aligned == THETA_STEPS - 1) begin
                           
                            accum_pipe_state <= PIPE_CALC_SUM_SHIFT;
                        end else begin
                           
                            accum_theta_idx  <= accum_theta_idx_aligned + 1; 
                            accum_pipe_state <= PIPE_CALC_SUM_SHIFT; 
                        end
                    end
                    PIPE_CALC_SUM_SHIFT: begin
                        hc_rho_sum_shifted_reg <= hc_rho_sum_shifted;
                        accum_theta_idx_p2     <= accum_theta_idx_p1; 
                        accum_pipe_state       <= PIPE_CALC_ADDR;
                    end
                    PIPE_CALC_ADDR: begin
                        accum_rd_addr_reg <= calculated_addr;
                        accum_pipe_state  <= PIPE_READ_HOUGH_ACCUM;
                    end
                    PIPE_READ_HOUGH_ACCUM: begin 
                        accum_read_data_reg <= hough_accum_dout_b;
                        accum_wr_addr_reg   <= accum_rd_addr_reg;
                        accum_pipe_state    <= PIPE_MODIFY_ACCUM;
                    end
                    PIPE_MODIFY_ACCUM: begin 
                        accum_pipe_state <= PIPE_WRITE_HOUGH_ACCUM;
                    end
                    PIPE_WRITE_HOUGH_ACCUM: begin 
                        accum_pipe_state <= PIPE_LOOP_THETA_ITER;
                    end
                    PIPE_LOOP_THETA_ITER: begin
                        if(accum_processing) begin 
                            if(accum_theta_idx_aligned == THETA_STEPS - 1) begin 
                                accum_processing <= 1'b0;
                                accum_pipe_state <= PIPE_FETCH_COORD; 
                            end else begin 
                                
                                accum_theta_idx_aligned <= accum_theta_idx;
                                accum_pipe_state <= PIPE_CALC_MUL_L0; 
                            end
                        end else begin
                            accum_pipe_state <= PIPE_FETCH_COORD;
                        end
                    end
                endcase
            end

            S_FIND_PEAKS: begin
                find_peaks_addr_reg <= find_theta_cnt * RHO_DIAMETER + find_rho_cnt;

               
                if ((find_theta_cnt_reg < THETA_DIVIDER) && (find_theta_cnt_reg >= THETA_LANE_MIN)) begin
                    if (hough_accum_dout_b_reg > peak1_val_nms) begin
                        peak1_val_nms     <= hough_accum_dout_b_reg;
                        peak1_theta_nms   <= find_theta_cnt_reg;
                        peak1_rho_idx_nms <= find_rho_cnt_reg;
                    end
                end

       
                if ((find_theta_cnt_reg >= THETA_DIVIDER) && (find_theta_cnt_reg <= THETA_LANE_MAX)) begin
                    if (hough_accum_dout_b_reg > peak2_val_nms) begin
                        peak2_val_nms     <= hough_accum_dout_b_reg;
                        peak2_theta_nms   <= find_theta_cnt_reg;
                        peak2_rho_idx_nms <= find_rho_cnt_reg;
                    end
                end

                if (find_rho_cnt == RHO_DIAMETER - 1) begin
                    find_rho_cnt <= 0;
                    if (find_theta_cnt == THETA_STEPS - 1) begin
                        peaks_done_flag <= 1'b1;
                        calc_step_cnt   <= 0;
                    end else begin
                        find_theta_cnt <= find_theta_cnt + 1;
                    end
                end else begin
                    find_rho_cnt <= find_rho_cnt + 1;
                end
            end

            S_CALCULATE_LINES: begin
                case(calc_step_cnt)
                    0: begin 
                        if (peak1_val_nms >= MIN_PEAK_THRESHOLD) begin
                            sorted_theta_right <= peak1_theta_nms;
                            sorted_rho_idx_right <= peak1_rho_idx_nms;
                        end else begin
                            sorted_theta_right <= 0;
                            sorted_rho_idx_right <= 0;
                        end

                        if (peak2_val_nms >= MIN_PEAK_THRESHOLD) begin
                            sorted_theta_left <= peak2_theta_nms;
                            sorted_rho_idx_left <= peak2_rho_idx_nms;
                        end else begin
                            sorted_theta_left <= 0;
                            sorted_rho_idx_left <= 0;
                        end
                        calc_step_cnt <= 1;
                    end
                    1: begin 
                        current_cand_theta_left    <= sorted_theta_left;
                        current_cand_rho_idx_left  <= sorted_rho_idx_left;
                        current_cand_theta_right   <= sorted_theta_right;
                        current_cand_rho_idx_right <= sorted_rho_idx_right;
                        calc_step_cnt <= 2;
                    end
                    2: begin 
                        if(!prev_final_left_valid_r) begin
                            final_theta_left <= current_cand_theta_left;
                            final_rho_idx_left <= current_cand_rho_idx_left;
                        end else if(is_left_lane_stable) begin
                            final_theta_left <= current_cand_theta_left;
                            final_rho_idx_left <= current_cand_rho_idx_left;
                        end else begin
                            final_theta_left <= prev_final_theta_left_r;
                            final_rho_idx_left <= prev_final_rho_idx_left_r;
                        end

                        if(!prev_final_right_valid_r) begin
                            final_theta_right <= current_cand_theta_right;
                            final_rho_idx_right <= current_cand_rho_idx_right;
                        end else if (is_right_lane_stable) begin
                            final_theta_right <= current_cand_theta_right;
                            final_rho_idx_right <= current_cand_rho_idx_right;
                        end else begin
                            final_theta_right <= prev_final_theta_right_r;
                            final_rho_idx_right <= prev_final_rho_idx_right_r;
                        end
                        calc_step_cnt <= 3;
                    end
                    3: begin 
                        wide_rho_true_left  <= $signed(final_rho_idx_left - RHO_MAX);
                        wide_rho_true_right <= $signed(final_rho_idx_right - RHO_MAX);
                        
                        calc_step_cnt <= 4;
                    end
                    4: begin 
                        wide_sec_lut_left_aligned     <= $signed(sec_val_bram_out_calc_left);
                        wide_neg_tan_lut_left_aligned <= $signed(neg_tan_val_bram_out_calc_left);
                        wide_sec_lut_right_aligned    <= $signed(sec_val_bram_out_calc_right);
                        wide_neg_tan_lut_right_aligned<= $signed(neg_tan_val_bram_out_calc_right);
                        calc_step_cnt <= 5;
                    end
                    5: begin 
                        left_term1  <= wide_rho_true_left  * wide_sec_lut_left_aligned;
                        left_term2  <= $signed(ROI_Y_START) * wide_neg_tan_lut_left_aligned;
                        right_term1 <= wide_rho_true_right * wide_sec_lut_right_aligned;
                        right_term2 <= $signed(ROI_Y_START) * wide_neg_tan_lut_right_aligned;
                        calc_step_cnt <= 6;
                    end
                    6: begin 
                        left_product  <= left_term1 + left_term2;
                        right_product <= right_term1 + right_term2;
                        calc_step_cnt <= 7;
                    end
                    7: begin
                        calc_done_flag <= 1'b1;
                        calc_step_cnt  <= 0;
                    end
                endcase
            end

            S_DRAW_FRAME: begin
                
                sin_val_draw_left_aligned <= $signed(sin_val_bram_out_draw_left);
                sin_val_draw_right_aligned <= $signed(sin_val_bram_out_draw_right);
                
                draw_y_sin_prod_left  <= $signed(y_out_cnt) * sin_val_draw_left_aligned; 
                draw_numerator_left   <= {wide_rho_true_left,  {COS_SIN_FRAC_BITS{1'b0}}} -
                                         {{1{draw_y_sin_prod_left[47]}}, draw_y_sin_prod_left};
                draw_product_left     <= draw_numerator_left  * wide_sec_lut_left_aligned; 
                x_lane_left           <= $signed(draw_product_left)  >>> TOTAL_DRAW_SHIFT;

                draw_y_sin_prod_right <= $signed(y_out_cnt) * sin_val_draw_right_aligned; 
                draw_numerator_right  <= {wide_rho_true_right, {COS_SIN_FRAC_BITS{1'b0}}} -
                                         {{1{draw_y_sin_prod_right[47]}}, draw_y_sin_prod_right};
                draw_product_right    <= draw_numerator_right * wide_sec_lut_right_aligned; 
                x_lane_right          <= $signed(draw_product_right) >>> TOTAL_DRAW_SHIFT;

                if (m_axis_tvalid && m_axis_tready) begin
                    if (x_out_cnt == IMG_WIDTH - 1) begin
                        x_out_cnt <= 0;
                        if (y_out_cnt == IMG_HEIGHT - 1) begin
                            y_out_cnt      <= 0;
                            draw_done_flag <= 1'b1;
                        end else begin
                            y_out_cnt <= y_out_cnt + 1;
                        end
                    end else begin
                        x_out_cnt <= x_out_cnt + 1;
                    end
                end
            end
        endcase
    end
end

endmodule
